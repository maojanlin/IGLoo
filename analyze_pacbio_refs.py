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
        if gene_name == "IGHD7-27":
            continue
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
    """ Gene pairs are (1, 0) instead of (0, 1) because it is not a pair in one segment but two consecutive segments """
    gene_0, min_dist_0 = find_closest_gene(pair_0, 1, list_gene_position, list_gene_name)
    gene_1, min_dist_1 = find_closest_gene(pair_1, 0, list_gene_position, list_gene_name)
    if abs(min_dist_0) < 50:
        if abs(min_dist_1) < 50:
            return True, gene_0, gene_1, min_dist_0, min_dist_1
        elif gene_1[3] == "V" and abs(min_dist_1) < 200:
            return True, gene_0, gene_1, min_dist_0, min_dist_1
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


def jump_and_fetch(fn_bam, list_alignment, seq_name):
    f_bam_sub = pysam.AlignmentFile(fn_bam)
    set_seg_info = set() # use set in case of segmental duplication making two segments start and the same position
    for c_id, SA_info in enumerate(list_alignment):
        for segment in f_bam_sub.fetch(SA_info[0], int(SA_info[1]), int(SA_info[1])+1):
            if seq_name == segment.query_name:
                query_seq = segment.query_alignment_sequence
                if query_seq != None:
                    foward_flag = segment.is_forward
                    start_pos = segment.reference_start # start position in genome coordiante
                    stop_pos  = segment.reference_end
                    cigar_tuples = segment.cigartuples
                    chr_name  =  SA_info[0]
                    set_seg_info.add((chr_name, start_pos, stop_pos, tuple(cigar_tuples), foward_flag, query_seq))
    return sorted(set_seg_info)


def seq_between_clips(cigar_tuples, seq_len, flag_forward):
    operation, length = cigar_tuples[0]
    if operation == 4 or operation == 5:
        pos_start = length
    else:
        pos_start = 0
    operation, length = cigar_tuples[-1]
    if operation == 4 or operation == 5:
        pos_end = seq_len - length
    else:
        pos_end = seq_len

    #print('----------------', len(cigar_tuples), flag_forward)
    if flag_forward:
        return pos_start, pos_end
    else:
        return seq_len-pos_start, seq_len-pos_end



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
    #print(list_target_idx[-1])
    #print(list_position)
    for idx in range(len(list_position) - 1):
        split_reads.append(sequence[list_position[idx]:list_position[idx+1]])
    return split_reads


def reverse_complement(seq):
    dict_reverse = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    r_seq = ""
    for ele in seq[::-1]:
        r_seq += dict_reverse[ele]
    return r_seq



def call_result_to_read_info(call_result, cigar_tuples, query_seq, gene_head, gene_tail, seg_st, seg_ed):
    """ call_result = call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name) """
    read_info = [[],[],[]]
    split_result = split_with_cigar([info[0] for info in call_result], cigar_tuples, query_seq)
    list_split_gene = []
    for split_info in call_result:
        list_split_gene.append(split_info[1])
        list_split_gene.append(split_info[2])
    # make positions
    list_split_gene = [gene_head] + list_split_gene + [gene_tail]
    list_seg_length = [len(split_seq) for split_seq in split_result]
    list_seg_pos    = []
    for idx in range(len(list_seg_length)):
        list_seg_pos.append(sum(list_seg_length[:idx+1]))
    list_seg_pos = [seg_st + seg_len for seg_len in list_seg_pos] if seg_st < seg_ed else [seg_st - seg_len for seg_len in list_seg_pos]
    list_seg_pos = [seg_st] + list_seg_pos
    for idx, split_seq in enumerate(split_result):
        read_info[0].append((list_split_gene[2*idx], list_split_gene[2*idx+1]))
        read_info[1].append((list_seg_pos[idx], list_seg_pos[idx+1]))
        read_info[2].append(split_seq)
    return read_info



def find_recombination(fn_bam, fn_bed, dict_read) -> dict:
    """
    dict_read_info = {}:
        - key: seq_name
        - value : ["type", [(gene0,gene1),(gene2,gene3),...], [(pos0,pos1),(pos2,pos3),...], [seq1,seq2,...]]
    
    type: "fir", "partial_fit", "non-fit", "Unrecombined", "D_read", "J_read"
    
    pos0, pos1 are the relative position on the whole read sequence
    """
    list_gene_position, list_gene_name, ref_name = read_bed(fn_bed)
    dict_gene_region = get_gene_region(list_gene_name, list_gene_position)

    J_region = dict_gene_region["J"]
    D_region = dict_gene_region["D"]
    V_region = dict_gene_region["V"]

    dict_read_info = {}
    
    f_bam = pysam.AlignmentFile(fn_bam, "rb")
    # Iterate from first J gene to first V gene
    for segment in f_bam.fetch(ref_name, J_region[0], V_region[0]):
        seq_name  = segment.query_name
        if dict_read_info.get(seq_name): # not to count supplemntary alignments twice
            continue
        dict_read_info[seq_name] = [] 

        start_pos = segment.reference_start # start position in genome coordiante
        stop_pos  = segment.reference_end
        cigar_tuples = segment.cigartuples
        flag_forward = segment.is_forward
        query_seq    = segment.query_alignment_sequence

        if not dict_read.get(seq_name):
            print("WARNING!", seq_name, "not in the fasta file")
            continue
        complete_seq = dict_read[seq_name] # segment.query_sequence doesn't include the hard clipped sequence
        if flag_forward:
            assert query_seq in complete_seq
        else:
            assert reverse_complement(query_seq) in complete_seq
        head_gene, min_dist_head = find_closest_gene(start_pos, 0, list_gene_position, list_gene_name)
        fst_gene,  min_dist_fst  = find_closest_gene(stop_pos,  1, list_gene_position, list_gene_name)

        # first check if there are supplementary alignment, if the split site is close to RSS, report, else defer for other reference genome
        if segment.has_tag("SA"):
            # make sure the "complete_seq" is sync with the pysam segments
            seg_st, seg_ed = seq_between_clips(cigar_tuples, len(complete_seq), flag_forward)
            if flag_forward:
                assert query_seq == complete_seq[seg_st:seg_ed]
            else:
                assert query_seq == reverse_complement( complete_seq[seg_ed:seg_st] )

            list_alignment = segment.get_tag("SA").split(";")[:-1]
            list_alignment = sorted([SA_tag.split(',') for SA_tag in list_alignment])
            list_seg_info  = jump_and_fetch(fn_bam, list_alignment, seq_name)
            
            # make sure the supplementary alignment is in sync with this reference genome, otherwise defer
            list_legit_seg  = []
            list_constant   = []
            list_non_fit    = []
            for seg_info in list_seg_info:
                contig_name = seg_info[0]
                seg_start   = seg_info[1]
                if contig_name == ref_name:
                    if seg_start >= stop_pos:
                        check_result = check_recomb_pair(stop_pos, seg_start, list_gene_position, list_gene_name)
                        if check_result[0]:
                            list_legit_seg.append(seg_info)
                        else:
                            list_non_fit.append(seg_info)
                    else: # segment before the J/D read
                        list_constant.append(seg_info)

            if list_legit_seg: # split and store to fo_s_fasta
                read_info = [[],[],[]]
                fst_head = head_gene if abs(min_dist_head) < 50 else None
                # check if there are internal deletions inside PRIMARY SEGMENTS
                call_result = call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name)
                if call_result:
                    read_info = call_result_to_read_info(call_result, cigar_tuples, query_seq, fst_head, fst_gene, seg_st, seg_ed)
                else: # normal aligned reads
                    read_info[0].append((fst_head, fst_gene))  # Add PRIMARY SEGMENT info
                    read_info[1].append((seg_st, seg_ed))      # Add PRIMARY SEGMENT info
                    read_info[2].append(None)      # Add PRIMARY SEGMENT info
                
                # first add the primary segment information
                for seg_info in list_legit_seg: # Fit split-segments
                    contig_name, seg_start, seg_stop, seg_tuples, seg_forward, seg_sequence = seg_info
                    seg_st, seg_ed = seq_between_clips(seg_info[3], len(complete_seq), seg_info[4])

                    gene_start, min_dist_start = find_closest_gene(seg_start, 0, list_gene_position, list_gene_name)
                    gene_stop, min_dist_stop   = find_closest_gene(seg_stop,  1, list_gene_position, list_gene_name)

                    info_start = gene_start if abs(min_dist_start) < 50 or (abs(min_dist_start) < 200 and gene_stop[3] == "V") else None
                    info_stop  = gene_stop  if abs(min_dist_stop)  < 50 else None
                    # check if there are internal split deletions in this segment
                    call_result = call_recomb_with_cigar(seg_tuples, seg_start, list_gene_position, list_gene_name)
                    if call_result:
                        read_split_info = call_result_to_read_info(call_result, seg_tuples, seg_sequence, info_start, info_stop, seg_st, seg_ed)
                        read_info[0] += read_split_info[0]
                        read_info[1] += read_split_info[1]
                        read_info[2] += read_split_info[2]
                    else:
                        read_info[0].append((info_start, info_stop))
                        read_info[1].append((seg_st, seg_ed))
                        read_info[2].append(None)

                for seg_info in list_constant: # Constant segments
                    contig_name, seg_start, seg_stop, seg_tuples, seg_forward, seg_sequence = seg_info
                    seg_st, seg_ed = seq_between_clips(seg_info[3], len(complete_seq), seg_info[4])
                    constant_gene, _ = find_closest_gene(seg_stop,  1, list_gene_position, list_gene_name)
                    
                    read_info[0].append((None, constant_gene))
                    read_info[1].append((seg_st, seg_ed))
                    read_info[2].append(None)
                for seg_info in list_non_fit: # non-fit segments:
                    contig_name, seg_start, seg_stop, seg_tuples, seg_forward, seg_sequence = seg_info
                    seg_st, seg_ed = seq_between_clips(seg_info[3], len(complete_seq), seg_info[4])
                    non_fit_gene, min_dist = find_closest_gene(seg_stop,  1, list_gene_position, list_gene_name)
                    
                    read_info[0].append((None, non_fit_gene)) if abs(min_dist) < 50 else read_info[0].append((None,None))
                    read_info[1].append((seg_st, seg_ed))
                    read_info[2].append(None)

                # assign to "partial fit" or "fit"
                if list_non_fit:
                    dict_read_info[seq_name] = ['partial_fit'] + read_info
                else:
                    dict_read_info[seq_name] = ['fit'] + read_info
                    
            else: # If there is no fitting spliting cites in the alignment
                dict_read_info[seq_name] = ["non-fit"]
        else:
            # if there is no supplementary alignment, check if there are >400 deletions similar to recombination
            call_result = call_recomb_with_cigar(cigar_tuples, start_pos, list_gene_position, list_gene_name)
            if call_result:
                seg_st, seg_ed = seq_between_clips(cigar_tuples, len(complete_seq), flag_forward)
                read_info = call_result_to_read_info(call_result, cigar_tuples, query_seq, None, None, seg_st, seg_ed)
                dict_read_info[seq_name] = ["fit"] + [read_info]
            else: # normal aligned reads
                last_element = [] if (len(complete_seq) - len(query_seq))<50 else [query_seq]
                if stop_pos > J_region[1] and start_pos < D_region[0]:
                    dict_read_info[seq_name] = ["Unrecombined",[],[],last_element]
                elif start_pos >= D_region[0]:
                    dict_read_info[seq_name] = ["D_read",[],[],last_element]
                else:
                    dict_read_info[seq_name] = ["J_read",[],[],last_element]

    return dict_read_info




        


def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-lbed', '--list_annotation_bed', required=True, nargs='+', help='the IGH annotations of the reference genomes.')
    parser.add_argument('-lbam', '--list_align_bam',      required=True, nargs='+', help='the alignment bam files of the IG reads to the reference genome.')
    parser.add_argument('-fasta', '--fasta',              required=True, help='the fasta file of the bams')
    parser.add_argument('-out',   '--output_report', help='the output report path')
    parser.add_argument('-out_fa', '--output_fasta', help='the output split read fasta')
    args = parser.parse_args(arguments)
    
    list_bed = args.list_annotation_bed
    list_bam = args.list_align_bam
    fn_fasta = args.fasta
    dict_read = read_fasta(fn_fasta)

    fn_out = args.output_report
    fn_output_fasta = args.output_fasta

    list_read_info = []
    for idx, fn_bam in enumerate(list_bam):
        fn_bed = list_bed[idx]
        list_read_info.append( find_recombination(fn_bam, fn_bed, dict_read) )

    for seq_name, sequence in dict_read.items():
        list_info = []
        for dict_read_info in list_read_info:
            if dict_read_info.get(seq_name):
                list_info.append(dict_read_info[seq_name])
        if list_info:
            category = list_info[0][0]
            for info in list_info[1:]:
                if info[0] != category:
                    print(seq_name)
                    print(list_info)



if __name__ == "__main__":
    main()


