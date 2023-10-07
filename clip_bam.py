import argparse
import pysam
import numpy as np
from analyze_recombine import read_bed, get_gene_region


def read_list(fn_list):
    set_target_read = set()
    fi = open(fn_list, 'r')
    for line in fi:
        set_target_read.add(line.strip())
    return set_target_read


def get_reverse_complement(input_seq):
    dict_reverse = {"A":"T", "T":"A", "C":"G", "G":"C", "a": "t", "t":"a", "c":"g", "g":"c"}
    rc_seq = ""
    for ele in input_seq[::-1]:
        rc_seq += dict_reverse[ele]
    return rc_seq


def cigar_read_len(cigar_string):
    bg_SH = 0
    ed_SH = 0
    read_len = 0
    current_len = 0
    bg_flag = True
    for ele in cigar_string:
        if ele.isdigit():
            current_len = int(ele) + current_len*10
        else:
            if ele == 'M' or ele == "I":
                read_len += current_len
            elif ele == 'S' or ele == "H":
                if bg_flag:
                    bg_SH = current_len
                else:
                    ed_SH = current_len
            current_len = 0
            bg_flag = False
    return bg_SH, read_len, ed_SH


def seperate_string_number(string):
    """borrow from Bryo Much and Nikaido"""
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i
        previous_character = i
        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups


def count_split_site(MD_tag, len_ref):
    list_diff = np.zeros(len_ref)
    list_idx = 0
    group_element = seperate_string_number(MD_tag)
    idx = 0
    while idx < len(group_element):
        ele = group_element[idx]
        if ele == '^': # deletion
            deletion = group_element[idx+1]
            idx += 1
            list_diff[list_idx:list_idx+len(deletion)] += 0.1
            list_idx += len(deletion)
        elif ele.isnumeric(): 
            list_idx += int(ele)
        else:
            list_diff[list_idx] += 1
            list_idx += 1
        idx += 1

    max_span = [0,0]
    max_len  = 0
    current_head = 0
    current_len  = 0
    first_block = sum(list_diff[:100])
    list_window = [first_block]
    for idx in range(1, len(list_diff)-100):
        list_window.append(first_block)
        if first_block < 11:
            if current_len == 0:
                current_head = idx
            current_len += 1
        else:
            if current_len > max_len:
                max_len = current_len
                max_span = [current_head, current_head + current_len]
            current_len = 0
        first_block += list_diff[idx+100]
        first_block -= list_diff[idx]
    if current_len > max_len:
        max_len = current_len
        max_span = [current_head, current_head + current_len]

    #print(max_span[0]+60, max_span[1]+40)
    if max_span[0] + 60 >= max_span[1]+40:
        print("WARNING! spanning too short!")
    return max_span[0]+60, max_span[1]+40
            






def clip_fasta(fn_bam, list_gene_position, list_gene_name, ref_name, \
               fn_fasta, fn_out, set_target_read):
    # standard read depleting
    fo = open(fn_out, 'w')
    fi = open(fn_fasta, 'r')

    
    # read clipping from bam information
    dict_gene_region = get_gene_region(list_gene_name, list_gene_position)
    J_region = dict_gene_region["J"]
    D_region = dict_gene_region["D"]
    
    """
    num_split = 0; num_no_split = 0
    f_bam = pysam.AlignmentFile(fn_bam)
    for segment in f_bam.fetch(ref_name, J_region[0]-12500, J_region[0]):
        seq_name  = segment.query_name
        if seq_name in set_target_read:
            continue
        if segment.has_tag("SA"):
            num_split += 1
        else:
            num_no_split +=1
    print(num_split, num_no_split)

    # discard constant splitting reads if there are non-spliting reads
    """
    set_constant_split = set()
    """
    if num_no_split > 5:
        f_bam = pysam.AlignmentFile(fn_bam)
        for segment in f_bam.fetch(ref_name, J_region[0]-12500, J_region[0]):
            seq_name  = segment.query_name
            if seq_name in set_target_read:
                continue
            if segment.has_tag("SA"):
                #print("###", seq_name, "split")
                query_seq    = segment.query_alignment_sequence
                set_constant_split.add(seq_name)
                fo.write('>' + seq_name + '/1\n')
                fo.write(query_seq  + '\n')
                
                original_seq = segment.query_sequence
                if len(original_seq) > len(query_seq):
                    fo.write('>' + seq_name + '/2\n')
                    fo.write(original_seq[original_seq.find(query_seq)+len(query_seq):] + '\n')
                #if flag_forward:
                #    print("###", seq_name, query_seq in original_seq, original_seq.find(query_seq), len(query_seq), len(original_seq))
                #else:
                #    print("RRR", seq_name, query_seq in original_seq, original_seq.find(query_seq), len(query_seq), len(original_seq))
            else:
                pass
                #print("###", seq_name, "NO")
            """
    
    ###### Adding the normal reads ######
    set_target_read.add("")
    dict_target_read = {}
    name = ""
    seq = ""
    for line in fi:
        if line[0] == '>':
            if name not in set_target_read and name not in set_constant_split:
                fo.write('>' + name + '\n')
                fo.write(seq + '\n')
            else:
                dict_target_read[name] = seq
            name = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
    # last sequence
    if name not in set_target_read and name not in set_constant_split:
        fo.write('>' + name + '\n')
        fo.write(seq + '\n')
    fi.close()
    ###### Adding the normal reads ######

    f_bam = pysam.AlignmentFile(fn_bam)
    for segment in f_bam.fetch(ref_name, J_region[0], J_region[1]):
        seq_name  = segment.query_name
        if seq_name not in set_target_read:
            continue
        start_pos = segment.reference_start # start position in genome coordiante
        stop_pos  = segment.reference_end
        cigar_tuples = segment.cigartuples
        cigar_string = segment.cigarstring
        cigar_states = segment.get_cigar_stats()
        flag_forward = not segment.is_reverse # forward: True, reverse: False
        MD_tag       = segment.get_tag("MD")
        query_seq    = segment.query_alignment_sequence
        #print(seq_name, segment.get_aligned_pairs()[-1], len(segment.query_alignment_sequence))

        """
        len_clip_L, len_clip_R = count_split_site(MD_tag, len(segment.get_reference_sequence()))
        len_clip_L += start_pos -1
        len_clip_R += start_pos -1
        flag_L = False
        split_L = 0
        split_R = len(query_seq) 
        for pair_info in segment.get_aligned_pairs():
            if len_clip_L == pair_info[1]:
                flag_L = True
            if len_clip_R == pair_info[1]:
                if pair_info[0] == None:
                    split_R = last_num
                else:
                    split_R = pair_info[0]
                break
            if pair_info[0] != None:
                if flag_L:
                    split_L = pair_info[0]
                    flag_L = False
                last_num = pair_info[0]

        fo.write(">" + seq_name + '/1\n')
        fo.write(query_seq[split_L:split_R] + '\n')
        """
        
        # list information of the supplementary alignments
        list_alignment = segment.get_tag("SA").split(";")[:-1]
        list_alignment = [ele.split(',') for ele in list_alignment]
        list_alignment = sorted(list_alignment, key=lambda x:int(x[1]))
        # Check if the read is in the fasta file
        if dict_target_read.get(seq_name) == None:
            print("Warning!", "'"+seq_name+"'", "not exist in fasta file.")
            continue
        
        """
        flag_strand_contradict = False
        if flag_forward:
            target_seq = dict_target_read[seq_name]
            for SA_info in list_alignment:
                if SA_info[2] != "+":
                    flag_strand_contradict = True
        else:
            target_seq = get_reverse_complement(dict_target_read[seq_name])
            for SA_info in list_alignment:
                if SA_info[2] != "-":
                    flag_strand_contradict = True
        if flag_strand_contradict:
            print("Flag strand contradict", seq_name)
            print(cigar_states)
            print(cigar_string)
            print(seq_name, start_pos, flag_forward, sum(cigar_states[0][:2]))
            for SA_info in list_alignment:
                print("=====", SA_info)
            continue"""


        
        
        seq_forward = dict_target_read[seq_name]
        seq_reverse = get_reverse_complement(seq_forward)
        flag_constant = False # check if there is constant gene recombination
        for SA_info in list_alignment:
            if SA_info[0] == ref_name and int(SA_info[1]) < start_pos: # it is a constant gene recombination
                flag_constant = True
                bg_SH, read_len, ed_SH = cigar_read_len(SA_info[3])
                
                fo.write(">" + seq_name + '/0\n')
                if SA_info[2] == "+": # forward
                    fo.write(seq_forward[bg_SH:bg_SH + read_len -300] + '\n')
                else:
                    fo.write(seq_reverse[bg_SH:bg_SH + read_len -300] + '\n')
        target_seq = seq_forward if flag_forward else seq_reverse
        bg_SH, read_len, ed_SH = cigar_read_len(cigar_string) # cutting from J sequence
        """
        if ed_SH != 0: # there are following sequence
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq[:-(ed_SH+60)] + '\n')
            if ed_SH > 300:
                fo.write(">" + seq_name + '/2\n')
                fo.write(target_seq[-(ed_SH-300):] + '\n')
        else:
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq + '\n')

        """
        if flag_constant: # if constant gene recombination, clip the beginning of the J gene sequence
            const_clip = 60
        else:
            const_clip = 0
            if bg_SH > 300:
                fo.write(">" + seq_name + '/0\n')
                fo.write(target_seq[:bg_SH] + '\n')
        if ed_SH != 0: # there are following sequence
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq[bg_SH+const_clip:-(ed_SH+60)] + '\n')
            if ed_SH > 300:
                fo.write(">" + seq_name + '/2\n')
                fo.write(target_seq[-(ed_SH-300):] + '\n')
        else:
            pass
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq[bg_SH+const_clip:] + '\n')

        """
        clip_dist = 0
        if list_alignment[0][0] == ref_name and int(list_alignment[0][1]) < start_pos:
            # Constant gene recombination here
            fo.write(">" + seq_name + '/0\n')
            fo.write(target_seq[:cigar_read_len(list_alignment[0][3])-300] + '\n')
            #print("\t", cigar_string)
            #print('\t', ','.join(list_alignment[0]), cigar_read_len(list_alignment[0][3]))
            assert(cigar_tuples[0][0] == 5 or cigar_tuples[0][0] == 4), cigar_string # make sure the sup-seq has clipping information
            j_gene_pos = cigar_tuples[0][1] # clipped length
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq[j_gene_pos+60:j_gene_pos+sum(cigar_states[0][:2])-clip_dist] + '\n')
        else:
            fo.write(">" + seq_name + '/1\n')
            fo.write(target_seq[:sum(cigar_states[0][:2])-clip_dist] + '\n')
        # right side of the split alignment
        if int(list_alignment[-1][1]) > D_region[1]:
            clip_dist = 300
        if cigar_read_len(list_alignment[-1][3]) > clip_dist:
            fo.write(">" + seq_name + '/2\n')
            fo.write(target_seq[-(cigar_read_len(list_alignment[-1][3])-clip_dist):] + '\n')
        """

    fo.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--annotation_bed', help='the IGH annotation of the reference genome.')
    parser.add_argument('-bam', '--align_bam', help='the alignment bam file of the IG reads to the reference genome.')
    parser.add_argument('-fasta', '--input_fasta', help='targeted original fasta containing unwanted reads.')
    parser.add_argument('-list', '--remove_list', help='list containing all the read name to be clipped.')
    parser.add_argument('-o', '--clipped_fasta', help='final fasta file without the unwanted reads.')
    args = parser.parse_args()

    fn_bed = args.annotation_bed
    fn_bam   = args.align_bam
    fn_fasta = args.input_fasta
    fn_list  = args.remove_list
    fn_out   = args.clipped_fasta
    
    list_gene_position, list_gene_name, ref_name = read_bed(fn_bed)

    set_target_read = read_list(fn_list)
    clip_fasta(fn_bam, list_gene_position, list_gene_name, ref_name, fn_fasta, fn_out, set_target_read)

