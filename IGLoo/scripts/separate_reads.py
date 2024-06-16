import argparse
import pysam

def read_fasta(fn_fasta):
    dict_read = {}
    name = ""
    seq  = ""
    with open(fn_fasta) as f:
        for line in f:
            if line[0] == '>':
                if name != "":
                    dict_read[name] = seq
                name = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
    dict_read[name] = seq
    return dict_read



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam',   '--input_bam', help='alignment file to H1 and H2 draft IGH locus.')
    parser.add_argument('-fasta', '--input_fasta', help='fasta file for the reads.')
    parser.add_argument('-out',   '--output_name', help='the output header of the H1 and H2 read fasta.')
    args = parser.parse_args(arguments)
    
    dict_read = read_fasta(args.input_fasta)
    
    dict_output = {'H1':{}, 'H2':{}}
    
    f_bam = pysam.AlignmentFile(args.input_bam, "rb")
    list_contigs = f_bam.references
    for contig_name in list_contigs:
        if 'H1' in contig_name:
            haplotype = 'H1'
        elif 'H2' in contig_name:
            haplotype = 'H2'
        else:
            continue

        for segment in f_bam.fetch(contig_name):
            q_name = segment.query_name
            if dict_read.get(q_name):
                dict_output[haplotype][q_name] = dict_read[q_name]

    output_header = args.output_name
    fo = open(output_header + '.H1.fa', 'w')
    for name, seq in dict_output['H1'].items():
        fo.write('>' + name + '\n')
        fo.write(seq + '\n')
    fo.close()
    fo = open(output_header + '.H2.fa', 'w')
    for name, seq in dict_output['H2'].items():
        fo.write('>' + name + '\n')
        fo.write(seq + '\n')
    fo.close()



    
if __name__ == "__main__":
    main()

