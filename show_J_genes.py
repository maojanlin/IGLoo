import argparse


def show_J_genes(fn_bed):
    f = open(fn_bed)
    dict_contig = {}
    for line in f:
        fields = line.strip().split()
        if len(fields) >= 3:
            contig = fields[0]
            gene = fields[3]
            if gene[:4] == "IGHJ":
                if dict_contig.get(contig):
                    dict_contig[contig].append(gene.split('*')[0])
                else:
                    dict_contig[contig] = [gene.split('*')[0]]
                    
    for contig, list_gene in dict_contig.items():
        for gene in list_gene:
            print(gene, end='\t')
        print("|||", end='\t')
    print()
    f.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--groub_bed_report', help='the group bed report file containing the grouping information of gAIRR-annotate')
    args = parser.parse_args()
    
    fn_bed = args.groub_bed_report
    show_J_genes(fn_bed)
