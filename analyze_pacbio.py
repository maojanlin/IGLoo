import argparse
import pysam
import subprocess
import os
import sys


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
    with open(fn_bed, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            fields = line.split()
            position = (int(fields[1]), int(fields[2]))
            list_gene_position.append(position)
            list_gene_name.append(fields[3].split('*')[0])
        chr_name = fields[0]
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


def get_gene_region(list_gene_name, list_gene_position):
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



def check_recomb_pair(pair_0, pair_1, list_gene_position, list_gene_name):
    gene_0, min_dist_0 = find_closest_gene(pair_0, 1, list_gene_position, list_gene_name)
    gene_1, min_dist_1 = find_closest_gene(pair_1, 0, list_gene_position, list_gene_name)
    if abs(min_dist_0) < 50 and abs(min_dist_1) < 50:
        return True, gene_0, gene_1, min_dist_0, min_dist_1
    else:
        return False, gene_0, gene_1, min_dist_0, min_dist_1
    

def call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name):
    list_recomb_candidate = []
    for idx, (operation, length) in enumerate(cigar_tuples):
        if operation == 2 and length > 400:
            list_recomb_candidate.append(idx)
    list_candidate_pair = []
    if list_recomb_candidate:
        position = start_pos
        for idx, (operation, length) in enumerate(cigar_tuples):
            if idx in list_recomb_candidate:
                list_candidate_pair.append((position, position+length))
            if operation == 4 or operation == 5 or operation == 1: # S or H or I
                continue
            elif operation == 0 or operation == 2: # M or D
                position += length
            else: # N or P or = or X or B
                print("WARNING!", seq_name, "with unsupported cigar string")
    list_recomb = []
    for idx, pair_info in enumerate(list_candidate_pair):
        check_result = check_recomb_pair(pair_info[0], pair_info[1], list_gene_position, list_gene_name)
        if check_result[0]:
            list_recomb.append([list_recomb_candidate[idx]] + list(check_result[1:]))
    return list_recomb


def read_fasta(fn_fasta):
    dict_read = {}
    seq = ""
    name = ""
    with open(fn_fasta) as f:
        for line in f:
            if line[0] == ">":
                if name != "":
                    dict_read[name] = seq
                name = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
        dict_read[name] = seq
    return dict_read



def jump_and_split(fn_bam, fo, list_alignment, seq_name, list_gene_position=None, list_gene_name=None):
    f_bam_sub = pysam.AlignmentFile(fn_bam)
    for c_id, SA_info in enumerate(list_alignment):
        for segment in f_bam_sub.fetch(SA_info[0], int(SA_info[1]), int(SA_info[1])+1):
            if seq_name == segment.query_name:
                query_seq = segment.query_alignment_sequence
                if query_seq != None:
                    if list_gene_position and list_gene_name:
                        cigar_tuples = segment.cigartuples
                        sequence = segment.query_alignment_sequence
                        call_result = call_recomb_with_cigar(cigar_tuples, segment.reference_start, \
                                                             list_gene_position, list_gene_name)
                        if call_result:
                            split_result = split_with_cigar([info[0] for info in call_result], cigar_tuples, sequence)
                            for idx, split_seq in enumerate(split_result):
                                fo.write('>' + seq_name + '/segment' + str(c_id+1) + '/sub' + str(idx) + SA_info[-1] + '_read\n')
                                fo.write(split_seq + '\n')
                        else:
                            fo.write('>' + seq_name + '/segment' + str(c_id+1) + '/' + SA_info[-1] + '_read\n')
                            fo.write(query_seq  + '\n')
                    else:
                        fo.write('>' + seq_name + '/segment' + str(c_id+1) + '/' + SA_info[-1] + '_read\n')
                        fo.write(query_seq  + '\n')
                break


def write_read(fo, seq_name, sequence):
    fo.write('>' + seq_name + '\n')
    fo.write(sequence + '\n')


def split_with_cigar(list_target_idx, cigar_tuples, sequence):
    list_position = []
    position = 0
    for idx in range(list_target_idx[-1]+1):
        (operation, length) = cigar_tuples[idx]
        if idx in list_target_idx:
            list_position.append(position)
        if operation == 4 or operation == 5 or operation == 2: # S or H or D
            continue
        elif operation == 0 or operation == 1: # M or I
            position += length
        else: # N or P or = or X or B
            position += length
            print("WARNING!", seq_name, "with unsupported cigar string")
    split_reads = []
    list_position = [0] + list_position + [-1]
    print(list_target_idx[-1])
    print(list_position)
    for idx in range(len(list_position) - 1):
        split_reads.append(sequence[list_position[idx]:list_position[idx+1]])
    return split_reads



def find_recombination(fn_bam, list_gene_name, list_gene_position, ref_name, dict_read, fn_out, fn_output_fasta, fn_defer_fasta):
    dict_gene_region = get_gene_region(list_gene_name, list_gene_position)
    if os.path.isfile(fn_output_fasta):
        fo_s_fasta = open(fn_output_fasta, 'a')
    else:
        fo_s_fasta = open(fn_output_fasta, 'w')
    fo_d_fasta = open(fn_defer_fasta, 'w')
    fo_rpt     = open(fn_out, 'w')

    J_region = dict_gene_region["J"]
    D_region = dict_gene_region["D"]
    V_region = dict_gene_region["V"]
    dict_VDJ_pairs = {tuple('Unrecombined'): []}
    dict_bad_reads = {}
    set_processed = set()
    
    f_bam = pysam.AlignmentFile(fn_bam, "rb")
    dict_record = {}
    for segment in f_bam.fetch(ref_name, J_region[0], V_region[0]):
        seq_name  = segment.query_name
        if seq_name in set_processed: # not to count supplemntary alignments twice
            continue
        set_processed.add(seq_name)

        start_pos = segment.reference_start # start position in genome coordiante
        stop_pos  = segment.reference_end
        cigar_tuples = segment.cigartuples
        fst_gene, min_dist_J = find_closest_gene(stop_pos, 1, list_gene_position, list_gene_name)

        # first check if there are supplementary alignment, if the split site is close to RSS, report, else defer for other reference genome
        if segment.has_tag("SA"):
            list_alignment = segment.get_tag("SA").split(";")[:-1]
            list_alignment = sorted([SA_info.split(',') for SA_info in list_alignment])
            list_const = []
            list_gene =  []
            list_SA_const = []
            list_SA_gene  = []
            print(seq_name)
            for SA_info in list_alignment:
                if SA_info[0] == ref_name:
                    SA_start = int(SA_info[1])
                    if SA_start < stop_pos:  # recombine between C and J
                        cigar_code = SA_info[3]
                        cigar_tuples = get_cigar_tuples(cigar_code)
                        SA_stop = get_stop_from_cigar(SA_start, cigar_tuples)
                        closest_C, min_dist_C = (find_closest_gene(SA_stop, 1, list_gene_position, list_gene_name))
                        print('\t', closest_C, min_dist_C)
                        list_const.append(closest_C)
                        list_SA_const.append(SA_info + ['C'])
                    else: # recombine between J/D/V
                        check_result = check_recomb_pair(stop_pos, SA_start, list_gene_position, list_gene_name)
                        if check_result[0]:
                            if list_gene == []:
                                list_gene.append(check_result[1])
                            list_gene.append(check_result[2])
                            print('-----', check_result)
                        else:
                            print('\t\t', check_result)
                        list_SA_gene.append(SA_info + [check_result[2][3]])
            if list_gene:
                call_result = call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name)
                if call_result:
                    split_result = split_with_cigar([info[0] for info in call_result], cigar_tuples, sequence)
                    for idx, split_seq in enumerate(split_result):
                        write_read(fo_s_fasta, seq_name + '/segment0/sub' + str(idx) + '/' + fst_gene[3] + '_read', split_seq)
                else:
                    write_read(fo_s_fasta, seq_name + '/segment0/' + fst_gene[3] + '_read', segment.query_alignment_sequence)
                jump_and_split(fn_bam, fo_s_fasta, list_SA_const, seq_name)
                jump_and_split(fn_bam, fo_s_fasta, list_SA_gene , seq_name, list_gene_position, list_gene_name)
            else: # defer read to fo_d_fasta
                write_read(fo_d_fasta, seq_name, dict_read[seq_name])
        else:
            # if there is no supplementary alignment, check if there are >1000 deletions similar to recombination
            call_result = call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name)
            if call_result:
                print(seq_name, call_result)
                split_result = split_with_cigar([info[0] for info in call_result], cigar_tuples, segment.query_alignment_sequence)
                for idx, split_seq in enumerate(split_result):
                    fo_s_fasta.write('>' + seq_name + '/segment' + str(idx) + '\n')
                    fo_s_fasta.write(split_seq + '\n')
            else: # normal aligned reads
                if stop_pos > J_region[1] and start_pos < D_region[0]:
                    print("Unrecombined read", seq_name)
                    write_read(fo_s_fasta, seq_name + '/Unrecombined', segment.query_alignment_sequence)
                elif start_pos >= D_region[0]:
                    print("Normal D read", seq_name)
                    write_read(fo_s_fasta, seq_name + '/D_read', segment.query_alignment_sequence)
                else:
                    print("Normal J read", seq_name)
                    write_read(fo_s_fasta, seq_name + '/J_read', segment.query_alignment_sequence)

    for seq_name, seq in dict_read.items():
        if seq_name not in set_processed:
            write_read(fo_s_fasta, seq_name, seq)

    fo_s_fasta.close()
    fo_d_fasta.close()





        


def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--annotation_bed',  required=True, help='the IGH annotation of the reference genome.')
    parser.add_argument('-bam', '--align_bam',       required=True, help='the alignment bam file of the IG reads to the reference genome.')
    parser.add_argument('-fasta', '--fasta',         help='the fasta file of the bam, generate one if not specified.')
    parser.add_argument('-out',  '--output_report',  required=True, help='the output report path')
    parser.add_argument('-out_fa',  '--output_fasta', help='the output split read fasta')
    parser.add_argument('-df_fa',  '--defered_fasta', help='the defered read fasta')
    args = parser.parse_args(arguments)
    
    fn_bed = args.annotation_bed
    fn_bam = args.align_bam
    fn_fasta = args.fasta
    if fn_fasta == None:
        fn_fasta = fn_bam + '.fa'
        subprocess.run(['samtools', 'fasta', fn_bam, '-o', fn_fasta, '-0', fn_fasta]) # '-0' for pacbio single reads
    dict_read = read_fasta(fn_fasta)

    fn_out = args.output_report
    fn_output_fasta = args.output_fasta
    fn_defer_fasta  = args.defered_fasta

    list_gene_position, list_gene_name, ref_name = read_bed(fn_bed)
    find_recombination(fn_bam, list_gene_name, list_gene_position, ref_name, dict_read, fn_out, fn_output_fasta, fn_defer_fasta)
    #find_overlap(fn_bam, list_gene_name, list_gene_position, ref_name, fn_out, flag_read_name, debug=True)



if __name__ == "__main__":
    main()

