import argparse
import pysam

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def find_junction(f_bam, f_fasta, target_spot):
    for segment in f_bam:
        flag = segment.flag
        if (flag & 4): # bitwise AND 4, segment unmapped
            continue
        # aligned read information
        ref_name     = segment.reference_name
        seq_name     = segment.query_name
        pos_start    = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        pos_end      = segment.reference_end
        cigar_tuples = segment.cigartuples
        mapq         = segment.mapping_quality
        read_seq     = segment.query_alignment_sequence # aligned sequence without SoftClip part

        if abs(pos_start - target_spot) < 10:
            # parameters 50, 100 are for the D, and J genes. Should be 300, 100 for V, J genes
            assert cigar_tuples[0][0] == 5 or cigar_tuples[0][0] == 4
            if (flag & 16): # if the sequence is reversed compliment
                read_length = f_fasta.get_reference_length(seq_name)
                spot_on_read = read_length - cigar_tuples[0][1]
                read_seq = f_fasta.fetch(reference=seq_name, start= spot_on_read - 50, end = spot_on_read + 100)
                print(">" + seq_name)
                print(read_seq)
            else:
                spot_on_read = cigar_tuples[0][1]
                read_seq = f_fasta.fetch(reference=seq_name, start= spot_on_read - 100, end = spot_on_read + 50)
                print(">" + seq_name)
                print(get_reverse_complement(read_seq))




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam',   '--alignment_file', help='the file containing the alignment and ,hopefully , junction information.')
    parser.add_argument('-fasta', '--read_file', help='the file containing the read sequence.')
    parser.add_argument('-spot', '--target_spot', type=int, help='the target splicing spot on the right side of the genome')
    args = parser.parse_args()
    
    fn_bam   = args.alignment_file
    fn_fasta = args.read_file
    target_spot = args.target_spot

    f_bam   = pysam.AlignmentFile(fn_bam)
    f_fasta = pysam.FastaFile(fn_fasta)

    find_junction(f_bam, f_fasta, target_spot)
