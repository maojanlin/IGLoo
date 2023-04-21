import argparse


list_IGHD = ['IGHD1-1','IGHD2-2','IGHD3-3','IGHD4-4','IGHD5-5','IGHD6-6','IGHD1-7','IGHD2-8','IGHD3-9','IGHD3-10','IGHD4-11','IGHD5-12','IGHD6-13','IGHD1-14','IGHD2-15','IGHD3-16', \
             'IGHD4-17','IGHD5-18','IGHD6-19','IGHD1-20','IGHD2-21','IGHD3-22','IGHD4-23','IGHD5-24','IGHD6-25','IGHD1-26','IGHD7-27']
list_IGHJ = ["IGHJ1P","IGHJ1","IGHJ2","IGHJ2P","IGHJ3","IGHJ4","IGHJ5","IGHJ3P","IGHJ6"]
list_IGLJ = ["IGLJ1","IGLJ2","IGLJ3","IGLJ4","IGLJ5","IGLJ6","IGLJ7"]
list_IGKJ = ["IGKJ1","IGKJ2","IGKJ3","IGKJ4","IGKJ5"]

def count_DJ(fn_bed,debug=True):
    dict_DJ = {}
    for gene in list_IGHD:
        dict_DJ[gene] = 0
    for gene in list_IGHJ:
        dict_DJ[gene] = 0
    for gene in list_IGLJ:
        dict_DJ[gene] = 0
    for gene in list_IGKJ:
        dict_DJ[gene] = 0

    f = open(fn_bed, 'r')
    dict_IG_contigs = {'IGH':{}, 'IGL': {}, 'IGK': {}}
    locus_flag = None
    for line in f:
        if line.strip() == 'IGH':
            locus_flag = 'IGH'
        elif line.strip() == 'IGL':
            locus_flag = 'IGL'
        elif line.strip() == 'IGK':
            locus_flag = 'IGK'

        if locus_flag == None:
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        contig_name = fields[0]
        allele_name = fields[3]
        gene_name   = allele_name.split('*')[0]
        
        if dict_DJ.get(gene_name) != None:
            dict_DJ[gene_name] += 1
    f.close()
    
    print(','.join(list_IGHD) + ',' + ','.join(list_IGHJ) + ',' + ','.join(list_IGLJ) + ',' + ','.join(list_IGKJ))
    for gene_name in list_IGHD:
        print(',', dict_DJ[gene_name], end="")
    for gene_name in list_IGHJ:
        print(',', dict_DJ[gene_name], end="")
    for gene_name in list_IGLJ:
        print(',', dict_DJ[gene_name], end="")
    for gene_name in list_IGKJ:
        print(',', dict_DJ[gene_name], end="")
    print()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--group_bed', help='the grouping bed file from gAIRR-annotate result.')
    args = parser.parse_args()
    
    fn_bed = args.group_bed

    print(fn_bed.split('/')[-2], ',', fn_bed.split('/')[-1].split('.')[-2], end="")
    count_DJ(fn_bed)
