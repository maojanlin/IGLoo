import argparse
import pysam
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
    read_len = 0
    current_len = 0
    for ele in cigar_string:
        if ele.isdigit():
            current_len = int(ele) + current_len*10
        else:
            if ele == 'M' or ele == "I":
                read_len += current_len
            current_len = 0
    return read_len


def clip_fasta(f_bam, list_gene_position, list_gene_name, ref_name, \
               fn_fasta, fn_out, set_target_read):
    # standard read depleting
    fo = open(fn_out, 'w')
    fi = open(fn_fasta, 'r')

    set_target_read.add("")
    dict_target_read = {}
    
    name = ""
    seq = ""
    for line in fi:
        if line[0] == '>':
            if name not in set_target_read:
                fo.write('>' + name + '\n')
                fo.write(seq + '\n')
            else:
                dict_target_read[name] = seq
            name = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
    # last sequence
    if name not in set_target_read:
        fo.write('>' + name + '\n')
        fo.write(seq + '\n')
    fi.close()

    # read clipping from bam information
    dict_gene_region = get_gene_region(list_gene_name, list_gene_position)
    J_region = dict_gene_region["J"]
    D_region = dict_gene_region["D"]
    print(D_region)
    
    for segment in f_bam.fetch(ref_name, J_region[0], J_region[1]):
        seq_name  = segment.query_name
        if seq_name not in set_target_read:
            continue
        start_pos = segment.reference_start # start position in genome coordiante
        stop_pos  = segment.reference_end
        cigar_tuples = segment.cigartuples
        cigar_string = segment.cigarstring
        cigar_states = segment.get_cigar_stats()
        flag_forward = not segment.is_reverse
        
        # list_alignment is like the applying
        list_alignment = segment.get_tag("SA").split(";")[:-1]
        list_alignment = [ele.split(',') for ele in list_alignment]
        list_alignment = sorted(list_alignment, key=lambda x:int(x[1]))
        #print(seq_name, start_pos, flag_forward, sum(cigar_states[0][:2]))
        #for SA_info in list_alignment:
        #    print("\t", ",".join(SA_info))
        flag_strand_contradict = False
        if dict_target_read.get(seq_name) == None:
            print("Warning!", "'"+seq_name+"'", "not exist in fasta file.")
            continue
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
            continue
        
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

    f_bam = pysam.AlignmentFile(fn_bam)
    set_target_read = read_list(fn_list)
    clip_fasta(f_bam, list_gene_position, list_gene_name, ref_name, fn_fasta, fn_out, set_target_read)

