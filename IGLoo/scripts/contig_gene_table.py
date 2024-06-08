import argparse
import numpy as np
import os



def solve_equivalent(list_gene_name, prev_gene, post_gene):
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
            append_element = [(1/2, "IGHV4-30-4"), (1/2, "IGHV4-31")]
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
    ###### For loop to read in the data ######
    for line in f:
        if line.strip() in ['IGH', 'IGL', 'IGK']:
            locus_flag = line.strip()
        if locus_flag == None:
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        
        contig_name = fields[0]
        allele_name = fields[3]
        position = int(fields[1])
        # Skipping the pseudogenes
        if ';' in allele_name: # dealing with identical allele name cases
            list_allele_name = allele_name.split(';')
            list_gene_name   = [name.split('*')[0] for name in list_allele_name]
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

        list_contig_allele.append((contig_name, allele_name, list_allele_name, list_gene_name, position))
    f.close()


    dict_contigs = {}
    dict_position = {}
    set_allele = set()
    set_gene   = set()
    ###### For loop to process the data ######
    for idx, info in enumerate(list_contig_allele):
        contig_name, allele_name, list_allele_name, list_gene_name, position = info
    
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
            append_element = solve_equivalent(list_gene_name, prev_gene, post_gene)
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
            dict_position[contig_name] += [(position, append_element[0][1])]
        else:
            dict_contigs[contig_name] = append_element
            dict_position[contig_name] = [(position, append_element[0][1])]
    
    return dict_contigs, set_allele, set_gene, dict_position


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


def report_break_points(dict_position, list_genes, fn_summary_bk, sample_name):
    """ Check V(D)J recombination and break points of the contigs """
    # make the report if the file isn't exist
    if os.path.isfile(fn_summary_bk):
        f = open(fn_summary_bk, 'a')
    else:
        f = open(fn_summary_bk, 'w')
        f.write('#sample_ID;VDJ_junctions;in_contig_SV;disjoint_bk;overlap_bk')

    # making the gene list table
    dict_list_genes = {}
    for idx, gene_name in enumerate(list_genes):
        dict_list_genes[gene_name] = idx
    set_VDJ_bk_point = set()
    
    # Finding all the VDJ junction genes on the contigs
    for contig_name, contig_info in dict_position.items():
        prev_pos  = contig_info[0][0]
        prev_gene = contig_info[0][1]
        for position, gene_name in contig_info:
            if prev_gene[3] != gene_name[3]:
                diff_gene = position - prev_pos
                if normal_VDJ_connection(prev_gene, gene_name) and (diff_gene > 10000 or "IGHV7-27" in (prev_gene, gene_name)):
                    pass
                else:
                    #print("\t\t", contig_name, prev_gene, gene_name, diff_gene)
                    #if prev_gene == "IGHV6-1" or gene_name == "IGHV6-1":
                    #    print("ATTENTION!!!!")
                    set_VDJ_bk_point.add(dict_list_genes[gene_name])
                    set_VDJ_bk_point.add(dict_list_genes[prev_gene])
            prev_pos  = position
            prev_gene = gene_name
    
    # Transform the genes into indecies
    list_index_genes = []
    for contig_name, contig_info in dict_position.items():
        list_index_genes.append([ dict_list_genes[info[1]] for info in contig_info ])
    # Finding all the remaining break points
    list_bk_points = []
    covered_genes  = np.zeros(len(list_genes)) # record the gene covering state of each contigs
    for indecies in list_index_genes:
        prev_index = indecies[0]
        st_index = 0
        for idx, index in enumerate(indecies):
            if abs(index - prev_index) > 5:
                if (prev_index in set_VDJ_bk_point) and (index in set_VDJ_bk_point): # the deletion due to VDJ recombination
                    #print("\tVDJ gap", index, prev_index)
                    pass
                else:
                    #print("\tPotential SV", index, prev_index)
                    list_bk_points.append(sorted((prev_index, index)))
                max_value = max(indecies[st_index:idx])
                min_value = min(indecies[st_index:idx])
                covered_genes[min_value:max_value+1] += 1
                st_index = idx
            prev_index = index
        max_value = max(indecies[st_index:idx+1])
        min_value = min(indecies[st_index:idx+1])
        covered_genes[min_value:max_value+1] += 1

    #print(covered_genes)
    #print("VDJ junction:", set_VDJ_bk_point)
    #print("inside contig SV:", list_bk_points)
    # Assigning each break points
    set_disjoint = set()
    set_overlap  = set()
    for indecies in list_index_genes:
        for index in [indecies[0], indecies[-1]]:
            if index not in [0, 94]:
                if covered_genes[index] == 0:
                    print("WARNING! uncovered gene in the list, there may be bug in the above code", index)
                elif covered_genes[index] == 1: 
                    set_disjoint.add(index)
                else:
                    set_overlap.add(index)
    #print("Disjoint break point:", set_disjoint)
    #print("Overlap  break point:", set_overlap)
    
    f.write('\n'+sample_name+';')
    f.write(','.join([list_genes[idx] for idx in sorted(set_VDJ_bk_point)]) + ';') # VDJ break points
    f.write(','.join([list_genes[ele[0]]+'&'+list_genes[ele[1]] for ele in list_bk_points]) + ';') # inside contig break points
    f.write(','.join([list_genes[idx] for idx in sorted(set_disjoint)]) + ';') # disjoint contig break points
    f.write(','.join([list_genes[idx] for idx in sorted(set_overlap)]))        # overlap  contig break points
    f.close()




def normal_VDJ_connection(gene_1, gene_2):
    """return True if it is normal, otherwise return False"""
    gene_pair = sorted([gene_1, gene_2])
    if gene_pair in [["IGHD1-26","IGHJ1P"], ["IGHD7-27","IGHJ1P"], ["IGHD7-27", "IGHJ1"], ["IGHD1-1", "IGHV6-1"]]:
        return True
    else:
        return False



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



def parse_gene(fn_target):
    list_genes=[]
    with open(fn_target) as f:
        for line in f:
            list_genes.append(line.strip())
    return list_genes





def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-target', '--target_genes', help='the file listing the genes for calling.', required=True)
    parser.add_argument('-bed1', '--group_bed_1', help='the first grouping bed file from gAIRR-annotate result.', required=True)
    parser.add_argument('-bed2', '--group_bed_2', help='the second grouping bed file from gAIRR-annotate result.')
    parser.add_argument('-out1', '--out_table_1', help='output table for the gene numbers in H1.')
    parser.add_argument('-out2', '--out_table_2', help='output table for the gene numbers in H2.')
    parser.add_argument('--summary_num', help='write the gene number summary result to the path')
    parser.add_argument('--summary_bk', help='write the break points summary result to the path')
    args = parser.parse_args(arguments)
    
    list_genes = parse_gene(args.target_genes)
    fn_bed_1 = args.group_bed_1
    fn_bed_2 = args.group_bed_2
    out_table_1 = args.out_table_1
    out_table_2 = args.out_table_2
    fn_summary_num = args.summary_num
    fn_summary_bk  = args.summary_bk

    dict_contigs_1, set_allele_1, set_gene_1, dict_position_1 = read_annotation(fn_bed_1, list_genes)
    contig_table(dict_contigs_1, list_genes, out_table_1)
    if fn_bed_2:
        dict_contigs_2, set_allele_2, set_gene_2, dict_position_2 = read_annotation(fn_bed_2, list_genes)
        contig_table(dict_contigs_2, list_genes, out_table_2)
        set_allele_1 = set_allele_1.union(set_allele_2)
        set_gene_1   = set_gene_1.union(set_gene_2)
    list_allele_num = count_VDJ(set_allele_1)
    list_gene_num   = count_VDJ(set_gene_1)

    # Summary of the individual
    if fn_summary_num:
        if os.path.isfile(fn_summary_num):
            f = open(fn_summary_num, 'a')
        else:
            f = open(fn_summary_num, 'w')
            f.write('#sample_ID,#V_gene,#D_gene,#J_gene,#V_allele,#D_allele,#J_allele')
        f.write('\n'+','.join([fn_bed_1.split('/')[-2]] + [str(num) for num in list_gene_num] + [str(num) for num in list_allele_num]))
        f.close()
    else:
        print(','.join([fn_bed_1.split('/')[-2]] + [str(num) for num in list_gene_num] + [str(num) for num in list_allele_num]))

    if fn_summary_bk:
        report_break_points(dict_position_1, list_genes, fn_summary_bk, fn_bed_1.split('/')[-2]+'.1')
        if fn_bed_2:
            report_break_points(dict_position_2, list_genes, fn_summary_bk, fn_bed_2.split('/')[-2]+'.2')
   


    

if __name__ == "__main__":
    main()
