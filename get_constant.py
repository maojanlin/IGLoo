import argparse
import pysam
import subprocess




def find_constant(f_bam, fn_fasta, fn_out):
    dict_gene = {}
    for segment in f_bam:
        gene_name = segment.query_name 
        flag = segment.flag
        if (flag & 4):
            print("WARNING!!", gene_name, "unmapped!")
            continue
        ref_name     = segment.reference_name
        pos_start    = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        pos_end      = segment.reference_end

        if dict_gene.get(ref_name):
            dict_gene[ref_name].append((gene_name, pos_start, pos_end))
        else:
            dict_gene[ref_name] = [(gene_name, pos_start, pos_end)]
    
    if len(dict_gene) > 1:
        print("WARNING!! There are", len(dict_gene), "segments contain constant gene!")
    if fn_out:
        f = open(fn_out, 'w')
        for ref_name, list_genes in dict_gene.items():
            min_pos = min([gene_pair[1] for gene_pair in list_genes])
            max_pos = max([gene_pair[2] for gene_pair in list_genes])
            #f.write('>' + list_genes[0][0] + '-' + list_gene[-1][0] + '-')
            #f.write(ref_name + ':' + str(min_pos) + '-' + str(max_pos) + '\n')
            constant_region = ref_name + ':' + str(min_pos) + '-' + str(max_pos)
            fasta_info = subprocess.run(['samtools', 'faidx', fn_fasta, constant_region], stdout=subprocess.PIPE)
            f.write(fasta_info.stdout.decode('utf-8'))
        f.close()
    else:
        for ref_name, list_genes in dict_gene.items():
            min_pos = min([gene_pair[1] for gene_pair in list_genes])
            max_pos = max([gene_pair[2] for gene_pair in list_genes])
            constant_region = ref_name + ':' + str(min_pos) + '-' + str(max_pos)
            subprocess.run(['samtools', 'faidx', fn_fasta, constant_region])
            #print(ref_name, min_pos, max_pos)
            #print(fasta_info)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam',   '--alignment_file', help='the file containing the alignment and ,hopefully , constant region information.')
    parser.add_argument('-fasta', '--ref_file', help='the file containing the assembly sequence.')
    parser.add_argument('-out',   '--output_file', help='the output file to contain the fetched constant region')
    args = parser.parse_args()
    
    fn_bam   = args.alignment_file
    fn_fasta = args.ref_file
    fn_out   = args.output_file

    f_bam   = pysam.AlignmentFile(fn_bam)

    find_constant(f_bam, fn_fasta, fn_out)
