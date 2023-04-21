import argparse

def check_IGH_exception(genes):
    set_chr15 = {"IGHV(II)-26-2", "IGHV(II)-1-1", "IGHV(III)-2-1"}
    set_chr16 = {"IGHV1-NL1", "IGHV(III)-16-1", "IGHV(II)-15-1", "IGHV(III)-13-1", 'IGHV(III)-11-1'}

    set_genes = set(genes)
    if set_genes.intersection(set_chr15) == set_genes:
        return "chr15_orphons"
    elif set_genes.intersection(set_chr16) == set_genes:
        return "chr16_orphons"
    else:
        return False


def count_IG_contigs(fn_bed, debug=False):
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
        
        if dict_IG_contigs[locus_flag].get(contig_name):
            dict_IG_contigs[locus_flag][contig_name].append(gene_name)
        else:
            dict_IG_contigs[locus_flag][contig_name] = [gene_name]
    f.close()

    num_IGH_contigs = 0
    for contig, genes in dict_IG_contigs['IGH'].items():
        flag_orphon = check_IGH_exception(genes) # skip if there is orphons
        if flag_orphon:
            if debug:
                print(flag_orphon)
            continue
        else:
            num_IGH_contigs += 1
    num_IGL_contigs = len(dict_IG_contigs['IGL'])
    num_IGK_contigs = len(dict_IG_contigs['IGK'])
    #print(' &', num_IGH_contigs, '&',  num_IGL_contigs,'&', num_IGK_contigs, '\\\\hline')
    print(',', num_IGH_contigs, ',',  num_IGL_contigs,',', num_IGK_contigs)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--group_bed', help='the grouping bed file from gAIRR-annotate result.')
    args = parser.parse_args()
    
    fn_bed = args.group_bed

    print(fn_bed.split('/')[-2], ',', fn_bed.split('/')[-1].split('.')[-2], end="")
    count_IG_contigs(fn_bed)
