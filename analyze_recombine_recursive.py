import argparse
import pysam
import subprocess


def print_combination(fn_bed):
    f = open(fn_bed, 'r')
    last_gene = ""
    last_name = ""
    for line in f:
        fields = line.split()
        if len(fields) < 4:
            last_gene = ""
            last_name = ""
            continue
        contig = fields[0]
        start  = int(fields[1])
        stop   = int(fields[2])
        gene_name = fields[3]
        current_gene = gene_name[3]

        if last_gene != "" and last_gene != current_gene:
            print(last_name, "---", gene_name)
        last_gene = current_gene
        last_name = gene_name

    f.close()


def read_bed(fn_bed):
    list_gene_position = []
    list_gene_name = []
    f = open(fn_bed, 'r')
    for line in f:
        if line[0] == '#':
            continue
        fields = line.split()
        position = (int(fields[1]), int(fields[2]))
        list_gene_position.append(position)
        list_gene_name.append(fields[3].split('*')[0])
    chr_name = fields[0]
    f.close()
    return list_gene_position, list_gene_name, chr_name


def get_cigar_tuples(cigar_string):
    cigar_tuples = []
    dict_cigar = {'M':0, 'I':1, 'D':2, 'S':4, 'H':5, '=':7, 'X':8}
    runs = 0
    for ele in cigar_string:
        if ele.isdigit():
            runs = runs*10 + int(ele)
        else:
            cigar_tuples.append((dict_cigar[ele], runs))
            runs = 0
    return cigar_tuples


def get_stop_from_cigar(
    read_start   :int,
    cigar_tuples :list
    ) -> int:
    """
    return the read stop position by read start and cigar
    """
    read_stop = read_start
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            read_stop += runs
        elif code == 1: # I
            pass
        elif code == 2: # D
            read_stop += runs
        elif code == 4 or code == 5: # S or H, pysam already parsed
            continue
        else:
            print("ERROR: unexpected cigar code in sequence", query_name)
    if code != 4 and code != 5:
        print("WARNING: not standard supplementary end!", cigar_tuples)
    return read_stop


def get_gene_region(list_gene_name, list_position):
    """
    each V, D, J return the maximum region expanded by the genes
    """
    dict_region = {"V":[3000000000,0], "D":[3000000000,0], "J":[3000000000,0], "C":[3000000000,0]}
    for idx, gene_name in enumerate(list_gene_name):
        gene_class = gene_name[3]
        if gene_class not in {"V", "D", "J"}:
            gene_class = "C"
        gene_region = list_gene_position[idx]
        if dict_region[gene_class][0] > gene_region[0]:
            dict_region[gene_class][0] = gene_region[0]
        if dict_region[gene_class][1] < gene_region[1]:
            dict_region[gene_class][1] = gene_region[1]
    return dict_region


def find_region(position, dict_gene_region):
    """
    return the region indicated by the gene type
    """
    if position > dict_gene_region["V"][1]:
        return("V+")
    elif position > dict_gene_region["V"][0]:
        return("V")
    elif position > dict_gene_region["D"][0]:
        return("D")
    elif position > dict_gene_region["J"][0]:
        return("J")
    elif position < dict_gene_region["C"][1]:
        return("C")
    else:
        return("C-J")


def find_closest_gene(position, flag_left, list_gene_position, list_gene_name):
    old_dist = 3000000000
    for idx, gene_region in enumerate(list_gene_position):
        current_dist = position - gene_region[flag_left]
        #print("..............................", position, current_dist, list_gene_name[idx], gene_region)
        if abs(current_dist) <= abs(old_dist):
            old_dist = current_dist
        else:
            return list_gene_name[idx-1], old_dist
    return list_gene_name[-1], old_dist


def find_JD_gap(dict_gene_region, list_gene_position, list_gene_name):
    JD_start = dict_gene_region["J"][1]
    JD_end   = -1
    for idx, gene_name in enumerate(list_gene_name):
        gene_class = gene_name[3]
        gene_region = list_gene_position[idx]
        if gene_class == "D" and gene_region[0] > JD_start:
            JD_end = gene_region[0]
            break
    return JD_start, JD_end


def find_overlap(fn_bam, list_gene_name, list_gene_position, ref_name, flag_read_name, debug=False):
    dict_gene_region = get_gene_region(list_gene_name, list_gene_position)

    J_region = dict_gene_region["J"]
    dict_VDJ_pairs = {tuple('Unrecombined'): []}
    
    f_alt_ref = open('./candidate_2_alt_ref.log', 'a')
    if debug:
        f = open('./count_distance.log', 'a')

    f_bam = pysam.AlignmentFile(fn_bam, "rb")
    for segment in f_bam.fetch(ref_name, J_region[0], J_region[1]):
        seq_name  = segment.query_name
        start_pos = segment.reference_start # start position in genome coordiante
        stop_pos  = segment.reference_end
        cigar_tuples = segment.cigartuples
        #print(seq_name)
        #print(find_closest_gene(stop_pos, 1, list_gene_position, list_gene_name))
        J_gene, min_dist_J = find_closest_gene(stop_pos, 1, list_gene_position, list_gene_name)
        if debug:
            f.write(fn_bam.split('.')[0] + "," + seq_name + ',' + J_gene + ',' + str(min_dist_J) + "\n")
        try:
            list_alignment = segment.get_tag("SA").split(";")[:-1]
            VD_gene = []
            C_gene = []
            for SA_tag in list_alignment:
                SA_info = SA_tag.split(',')
                if SA_info[0] != ref_name:
                    continue
                else:
                    SA_start = int(SA_info[1])
                    if find_region(SA_start, dict_gene_region) in ("D", "V"):
                        #print(find_region(SA_start, dict_gene_region), find_closest_gene(SA_start, 1, list_gene_position, list_gene_name))
                        # TODO can still be unrecombined sequences
                        # TODO examine for RSS sequences
                        closest_VD, min_dist_VD = (find_closest_gene(SA_start, 0, list_gene_position, list_gene_name))
                        VD_gene.append(closest_VD)
                        # TODO maybe add an option
                        if (min_dist_J > 0) or (min_dist_VD < 0) or abs(min_dist_J) < -500 or abs(min_dist_VD) > 500:
                            print(seq_name, "-----------------------------", J_gene, min_dist_J)
                            print(seq_name, '+++++++++++++++++++++++++++++', closest_VD, min_dist_VD)
                        if (min_dist_VD < -100) or (min_dist_VD) > 100:
                            f_alt_ref.write(fn_bam.split('.')[0] + "," + seq_name + ',' + closest_VD + ',' + str(min_dist_VD) + "\n")
                        if debug:
                            f.write(fn_bam.split('.')[0] + "," + seq_name + ',' + closest_VD + ',' + str(min_dist_VD) + "\n")
                    elif find_region(SA_start, dict_gene_region) in ("C"):
                        cigar_code = SA_info[3]
                        cigar_tuples = get_cigar_tuples(cigar_code)
                        SA_stop = get_stop_from_cigar(SA_start, cigar_tuples)
                        closest_VD, min_dist_VD = (find_closest_gene(SA_stop, 1, list_gene_position, list_gene_name))
                        if closest_VD[3] == "J":
                            continue
                        C_gene.append(closest_VD)
                        # TODO maybe add an option
                        if (min_dist_J > 0) or (min_dist_VD < 0) or abs(min_dist_J) < -500 or abs(min_dist_VD) > 500:
                            print(seq_name, "---------", J_gene, min_dist_J)
                            print(seq_name, '+++++++++', closest_VD, min_dist_VD)
                        #if (min_dist_VD < 0) or (min_dist_VD) > 4000:
                        #    f_alt_ref.write(fn_bam.split('.')[0] + "," + seq_name + ',' + closest_VD + ',' + str(min_dist_VD) + "\n")
                        if debug:
                            f.write(fn_bam.split('.')[0] + "," + seq_name + ',' + closest_VD + ',' + str(min_dist_VD) + "\n")
            if len(VD_gene) > 0:
                if (min_dist_J > 100) or (min_dist_J) < -100:
                    f_alt_ref.write(fn_bam.split('.')[0] + "," + seq_name + ',' + J_gene + ',' + str(min_dist_J) + "\n")
                if dict_VDJ_pairs.get(tuple(C_gene + [J_gene] + VD_gene)):
                    dict_VDJ_pairs[tuple(C_gene + [J_gene] + VD_gene)].append(seq_name)
                else:
                    dict_VDJ_pairs[tuple(C_gene + [J_gene] + VD_gene)] = [seq_name]
            elif stop_pos > J_region[1]: # The supplementary segments are all before J genes
                dict_VDJ_pairs[tuple('Unrecombined')].append(seq_name)
        except:
            #print("No SA")
            if stop_pos > J_region[1]:
                dict_VDJ_pairs[tuple('Unrecombined')].append(seq_name)
    for combination, list_seq_name in sorted(dict_VDJ_pairs.items()):
        if combination == tuple('Unrecombined'):
            print('Unrecombined', len(list_seq_name), sep='\t')
        else:
            print('---'.join(combination), len(list_seq_name), sep='\t')
        if flag_read_name:
            for read_name in list_seq_name:
                print(read_name)
    
    f_alt_ref.close()
    if debug:
        f.close()
    
    # show the average read depth between the J-D gap
    # get the J-D gap position
    JD_start, JD_end = find_JD_gap(dict_gene_region, list_gene_position, list_gene_name)
    gap_region = "chr14:" + str(JD_start) + '-' + str(JD_end)
    # call samtools to get the read depth
    read_depth = subprocess.run(['samtools', 'depth', '-r', gap_region, fn_bam], stdout=subprocess.PIPE)
    sum_read_depth = 0
    for line in read_depth.stdout.decode('utf-8').split('\n')[:-1]:
        sum_read_depth += int(line.split()[2])
    print("Average Read Depth between J-D gap:", sum_read_depth/ (JD_end - JD_start))







if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-lb', '--list_bed', nargs='+', required=True, help='the list of the IGH annotations of the reference genomes')
    parser.add_argument('-la', '--list_bam', nargs='+', required=True, help='the list of the alignment bam files')
    parser.add_argument('-bed', '--annotation_bed', help='the IGH annotation of the reference genome.')
    parser.add_argument('-bam', '--align_bam',      help='the alignment bam file of the IG reads to the reference genome.')
    parser.add_argument('-name', '--print_name_flag', action='store_true', help='print the seq name supporting a recombination type.')
    args = parser.parse_args()
    
    fn_bed = args.annotation_bed
    fn_bam = args.align_bam
    flag_read_name = args.print_name_flag

    list_gene_position, list_gene_name, ref_name = read_bed(fn_bed)
    find_overlap(fn_bam, list_gene_name, list_gene_position, ref_name, flag_read_name, debug=False)

