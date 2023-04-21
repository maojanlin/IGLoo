import argparse


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




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--group_bed', help='the grouping bed file from gAIRR-annotate result.')
    args = parser.parse_args()
    
    fn_bed = args.group_bed

    #print(fn_bed.split('/')[-2], ',', fn_bed.split('/')[-1].split('.')[-2], end="")
    print_combination(fn_bed)
