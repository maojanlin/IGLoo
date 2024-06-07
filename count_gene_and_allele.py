import argparse
from contig_gene_table import read_annotation, parse_gene


def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-target', '--target_genes', help='the file listing the genes for calling.')
    parser.add_argument('-bed1', '--group_bed_1', help='the first grouping bed file from gAIRR-annotate result.')
    parser.add_argument('-bed2', '--group_bed_2', help='the second grouping bed file from gAIRR-annotate result.')
    args = parser.parse_args(arguments)
    
    fn_bed_1 = args.group_bed_1
    fn_bed_2 = args.group_bed_2
    list_genes = parse_gene(args.target_genes)

    dict_contigs_1, set_allele_1, set_gene_1 = read_annotation(fn_bed_1)
    dict_contigs_2, set_allele_2, set_gene_2 = read_annotation(fn_bed_2)

    total_allele_num = len(set_allele_1.union(set_allele_2))
    total_gene_num = len(set_gene_1.union(set_gene_2))

    print(','.join([fn_bed_1.split('/')[-2], str(total_gene_num), str(total_allele_num)]))


if __name__ == "__main__":
    main()
