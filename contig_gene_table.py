import argparse
import os


def solve_duplicate(list_gene_name, prev_gene, post_gene):
    set_neighbor = set(prev_gene).union(set(post_gene))
    if list_gene_name == ["IGHD5-18","IGHD5-5"]: #IGHD5-18, IGHD5-5
        if "IGHD6-6" in set_neighbor or 'IGHD4-4' in set_neighbor:
            append_element = [(1, "IGHD5-5")]
        elif "IGHD6-19" in set_neighbor or 'IGHD4-17' in set_neighbor:
            append_element = [(1, "IGHD5-18")]
        else:
            append_element = [(1/2, "IGHD5-18"), (1/2, "IGHD5-5")]
    elif list_gene_name == ["IGHV4-30-4","IGHV4-31"]: # IGHV4-30-4, IGHV4-31
        if "IGHV3-33" in set_neighbor and "IGHV3-30-5" in set_neighbor:
            append_element = [(1, "IGHV4-31")]
        elif "IGHV3-33" in set_neighbor and "IGHV3-30-3" in set_neighbor:
            append_element = [(1, "IGHV4-30-4")]
        elif "IGHV3-30-5" in set_neighbor and "IGHV3-30-3" in set_neighbor:
            append_element = [(1, "IGHV4-30-4")]
        else:
            append_element = [(1/2, "IGHD4-30-4"), (1/2, "IGHD4-31")]
    elif list_gene_name == ["IGHV3-30", "IGHV3-30-5"]: # IGHV3-30, IGHV3-30-5
        if "IGHV4-28" in set_neighbor:
            append_element = [(1, "IGHV3-30")]
        else:
            append_element = [(1, "IGHV3-30-5")]
    elif list_gene_name == ["IGHV3-30", "IGHV3-30-3"]: # IGHV3-30, IGHV3-30-3
        if "IGHV4-28" in set_neighbor:
            append_element = [(1, "IGHV3-30")]
        else:
            append_element = [(1, "IGHV3-30-3")]
    else:
        append_element = [(1/len(list_gene_name), name) for name in list_gene_name]
    return append_element



def read_annotation(fn_bed, list_genes):
    f = open(fn_bed, 'r')
    locus_flag = None
    list_contig_allele = []
    ###### For loop to take in the data ######
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
        if ';' in allele_name:
            list_allele_name = allele_name.split(';')
            list_gene_name   = [name.split('*')[0] for name in list_allele_name]
            # skip the pseudogenes
            flag_functional = False
            for gene_name in list_gene_name:
                if gene_name in list_genes:
                    flag_functional = True
            if flag_functional == False:
                continue
        else:
            gene_name = allele_name.split('*')[0]
            if gene_name not in list_genes:
                continue
            list_allele_name = [allele_name]
            list_gene_name   = [gene_name]

        list_contig_allele.append((contig_name, allele_name, list_allele_name, list_gene_name))
    f.close()


    dict_contigs = {}
    set_allele = set()
    set_gene   = set()
    ###### For loop to process the data ######
    for idx, info in enumerate(list_contig_allele):
        contig_name, allele_name, list_allele_name, list_gene_name = info
    
        if len(list_gene_name) > 1:
            if idx == 0:
                prev_gene = ["None"]
            else:
                prev_gene = list_contig_allele[idx-1][3]
            if idx == len(list_contig_allele) -1:
                post_gene = ["None"]
            else:
                post_gene = list_contig_allele[idx+1][3]
            
            # only solve the identical gene IGHD5-18 and IGHD5-5,
            # the duplicate genes IGHV1-69, IGHV2-70, and IGHV3-23 are resolved later in function contig_table()
            append_element = solve_duplicate(list_gene_name, prev_gene, post_gene)
            for gene_name in [ele[1] for ele in append_element]:
                set_gene.add(gene_name)
                for allele_name in list_allele_name:
                    if allele_name.split('*')[0] == gene_name:
                        set_allele.add(allele_name)
        else:
            gene_name = list_gene_name[0]
            append_element = [(1, gene_name)]
            set_allele.add(allele_name)
            set_gene.add(gene_name)
        
        if dict_contigs.get(contig_name):
            dict_contigs[contig_name] += append_element
        else:
            dict_contigs[contig_name] = append_element
    
    return dict_contigs, set_allele, set_gene



def parse_gene(fn_target):
    list_genes=[]
    f = open(fn_target)
    for line in f:
        list_genes.append(line.strip())
    return list_genes



def curate_duplicate(dict_genes):
    #dict_duplicate = {"IGHV1-68", "IGHV1-69", "IGHV2-70", "IGHV3-23"}
    dict_duplicate = {"IGHV1-69":0, "IGHV2-70":0, "IGHV3-23":0}
    for dup_gene in dict_duplicate.keys():
        number = dict_genes[dup_gene] + dict_genes[dup_gene+"D"]
        dict_duplicate[dup_gene] = number

    for dup_gene, number in dict_duplicate.items():
        if number == 1:
            dict_genes[dup_gene] = 1
            dict_genes[dup_gene+"D"] = 0
        elif number >= 2:
            dict_genes[dup_gene] = number/2
            dict_genes[dup_gene+"D"] = number/2

    return dict_duplicate



def contig_table(dict_contigs, list_genes, out_table):
    """
    dict_genes: gene_count dict for each contig
    """
    if out_table:
        fo = open(out_table, 'w')
        fo.write(','.join(["contig"] + list_genes))
    else:
        print(','.join(["contig"] + list_genes))
    for contig_name, contig_genes in dict_contigs.items():
        dict_genes = {gene:0 for gene in list_genes}
        for gene_info in contig_genes:
            if dict_genes.get(gene_info[1]) != None:
                dict_genes[gene_info[1]] += gene_info[0]
        curate_duplicate(dict_genes) # curate the duplicate genes IGHV1-69, IGHV2-70, and IGHV3-23
        list_gene_number = [dict_genes[ele] for ele in list_genes]
        if sum(list_gene_number) >= 1:
            if out_table:
                fo.write('\n' + ','.join([contig_name] + [str(num) for num in list_gene_number]))
            else:
                print(','.join([contig_name] + [str(num) for num in list_gene_number]))
    
    if out_table:
        fo.close()


def count_VDJ(gene_names):
    num_V, num_D, num_J = 0, 0, 0
    for name in gene_names:
        if name[3] == "V":
            num_V += 1
        elif name[3] == "D":
            num_D += 1
        elif name[3] == "J":
            num_J += 1
        else:
            print("WARNING, unconventional gene/allele name:", name)
    return num_V, num_D, num_J






def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-target', '--target_genes', help='the file listing the genes for calling.', required=True)
    parser.add_argument('-bed1', '--group_bed_1', help='the first grouping bed file from gAIRR-annotate result.', required=True)
    parser.add_argument('-bed2', '--group_bed_2', help='the second grouping bed file from gAIRR-annotate result.')
    parser.add_argument('-out1', '--out_table_1', help='output table for the gene numbers in H1.')
    parser.add_argument('-out2', '--out_table_2', help='output table for the gene numbers in H2.')
    parser.add_argument('--summary', help='write the summary result to the path')
    args = parser.parse_args(arguments)
    
    list_genes = parse_gene(args.target_genes)
    fn_bed_1 = args.group_bed_1
    fn_bed_2 = args.group_bed_2
    out_table_1 = args.out_table_1
    out_table_2 = args.out_table_2
    fn_summary  = args.summary

    dict_contigs_1, set_allele_1, set_gene_1 = read_annotation(fn_bed_1, list_genes)
    contig_table(dict_contigs_1, list_genes, out_table_1)
    if fn_bed_2:
        dict_contigs_2, set_allele_2, set_gene_2 = read_annotation(fn_bed_2, list_genes)
        contig_table(dict_contigs_2, list_genes, out_table_2)
        set_allele_1 = set_allele_1.union(set_allele_2)
        set_gene_1   = set_gene_1.union(set_gene_2)
    list_allele_num = count_VDJ(set_allele_1)
    list_gene_num   = count_VDJ(set_gene_1)

    #print(','.join([fn_bed.split('/')[-2], fn_bed.split('/')[-1].split('.')[-2], str(len(set_gene)), str(len(set_allele))]))
    if fn_summary:
        if os.path.isfile(fn_summary):
            f = open(fn_summary, 'a')
        else:
            f = open(fn_summary, 'w')
            f.write('#sample_ID,#V_gene,#D_gene,#J_gene,#V_allele,#D_allele,#J_allele')
        f.write('\n'+','.join([fn_bed_1.split('/')[-2]] + [str(num) for num in list_gene_num] + [str(num) for num in list_allele_num]))
        f.close()
    else:
        print(','.join([fn_bed_1.split('/')[-2]] + [str(num) for num in list_gene_num] + [str(num) for num in list_allele_num]))
   


    

if __name__ == "__main__":
    main()
