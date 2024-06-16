import subprocess
import argparse
import numpy as np
import math



def parse_csv(csv_file):
    dict_contigs_genes = {}
    with open(csv_file) as f:
        header = f.readline()
        header = header.strip().split(',')[1:]
        for line in f:
            fields = line.strip().split(',')
            dict_contigs_genes[fields[0]] = [float(ele) for ele in fields[1:]]

    return header, dict_contigs_genes


def get_range(list_genes):
    bg = -1
    ed = -1
    for idx, ele in enumerate(list_genes):
        if ele != 0:
            if bg == -1:
                bg = idx
            ed = idx
    return bg, ed


def sv_range(header):
    """hard code the deletion and complex SV range"""
    deletion = [[] for idx in range(6)]
    deletion[0] = ['IGHD2-8', 'IGHD1-7', 'IGHD6-6', 'IGHD5-5', 'IGHD4-4', 'IGHD3-3']
    deletion[1] = ['IGHV7-4-1']
    deletion[2] = ['IGHV3-23D']
    deletion[3] = ['IGHV4-30-2', 'IGHV3-30-3', 'IGHV4-30-4', 'IGHV3-30-5', 'IGHV4-31']
    deletion[4] = ['IGHV4-38-2', 'IGHV3-43D', 'IGHV3-38-3', 'IGHV1-38-4']
    deletion[5] = ['IGHV2-70D', 'IGHV1-69-2', 'IGHV1-69D']
    complex_sv = [[],[]]
    complex_sv[0]  = ['IGHV1-8', 'IGHV3-9']
    complex_sv[1]  = ['IGHV3-64D', 'IGHV5-10-1']

    deletion_range = [[] for idx in range(len(deletion))]
    for idx, list_gene in enumerate(deletion):
        list_range = []
        for idy, gene_name in enumerate(header):
            if gene_name in list_gene:
                list_range.append(idy)
        deletion_range[idx] = [min(list_range), max(list_range)]
    
    complex_range = [[],[]]
    for idx in (0, 1):
        list_range = []
        for idy, gene_name in enumerate(header):
            if gene_name in complex_sv[idx]:
                list_range.append(idy)
        complex_range[idx] = [min(list_range), max(list_range)]
    return deletion_range, complex_range


def report_haplotype(list_genes, deletion_range, complex_range):
    """haplotype unknown: -1, 1copy: 1, 2copy: 2. complex in a similar fashion"""
    bg, ed = get_range(list_genes) 
    list_haplotype = []
    for sv_bg, sv_ed in deletion_range:
        if bg < sv_bg and ed > sv_ed:
            haplotype = sum(list_genes[sv_bg:sv_ed+1]) / (sv_ed-sv_bg+1)
            list_haplotype.append(round(haplotype))
            #list_haplotype.append(math.ceil(haplotype))
        else:
            list_haplotype.append(-1)

    complex_0 = sum(list_genes[complex_range[0][0]:complex_range[0][1]+1])
    complex_1 = sum(list_genes[complex_range[1][0]:complex_range[1][1]+1])
    if complex_0 > 0 and complex_1 == 0:
        list_haplotype.append(0)
    elif complex_1 > 0 and complex_0 == 0:
        list_haplotype.append(1)
    elif complex_0 > 0 and complex_1 > 0:
        print("WARNING! complex SV got both haplotypes in one contig!")
        list_haplotype.append(-1)
    else:
        list_haplotype.append(-1)
    return list_haplotype


def check_consistent(list_haplotype):
    num_contigs = len(list_haplotype)
    list_consensus = []
    flag_consist = True
    for idx, sv in enumerate(list_haplotype[0]):
        hap = -1
        list_consensus.append(hap)
        for idy in range(num_contigs):
            current_hap = list_haplotype[idy][idx]
            if current_hap == -1:
                continue
            elif hap == -1:
                hap = current_hap
            else:
                if current_hap == hap:
                    continue
                else:
                    flag_consist = False
                    hap = -1
                    break
        list_consensus[-1] = hap
    return flag_consist, list_consensus


def report_outlier(array_haplotype):
    list_diff = []
    for idx in range(len(array_haplotype[0])):
        sv_hap = array_haplotype[:,idx]
        target = sv_hap[sv_hap != -1]
        if len(target) <= 1:
            continue
        elif len(set(target)) == 1:
            continue
        else:
            outlier_hap = set()
            for idx, ele in enumerate(np.bincount(target)):
                if ele == 1:
                    outlier_hap.add(idx)
            for idx, ele in enumerate(sv_hap):
                if ele in outlier_hap:
                    list_diff.append(idx)
    return set(list_diff)


def sort_candidate_outlier(dict_contigs_genes, list_diff):
    sorted_contigs_genes = sorted(dict_contigs_genes.items())
    list_len_idx_name = []
    for idx in list_diff:
        name, list_genes = sorted_contigs_genes[idx]
        bg, ed = get_range(list_genes)
        list_len_idx_name.append((ed-bg, idx, name))
    return sorted(list_len_idx_name)


def reconcile(list_csv1_haplotype, list_csv2_haplotype, dict_contigs_genes_1, dict_contigs_genes_2):
    array_haplotype_1 = np.array(list_csv1_haplotype)
    array_haplotype_2 = np.array(list_csv2_haplotype)

    list_diff_1 = report_outlier(array_haplotype_1)
    list_diff_2 = report_outlier(array_haplotype_2)
    #print(array_haplotype_1)
    #print(array_haplotype_2)
    #print("........................")
    #print(list_diff_1)
    #print(list_diff_2)
    #print("........................")

    list_len_idx_name_1 = sort_candidate_outlier(dict_contigs_genes_1, list_diff_1)
    list_len_idx_name_2 = sort_candidate_outlier(dict_contigs_genes_2, list_diff_2)
    if not list_len_idx_name_1:
        list_len_idx_name_1 = [(-1,-1,-1)]
    if not list_len_idx_name_2:
        list_len_idx_name_2 = [(-1,-1,-1)]
    #print(list_len_idx_name_1)
    #print(list_len_idx_name_2)

    for _, idx_1, name_1 in list_len_idx_name_1:
        if idx_1 == -1:
            array_outlier_1 = (np.ones(len(array_haplotype_1[0]))*-1).astype(int)
            array_remain_1  = array_haplotype_1
        else:
            array_outlier_1 = array_haplotype_1[idx_1]
            array_remain_1  = np.vstack([array_haplotype_1[:idx_1,:], array_haplotype_1[idx_1+1:,:]])

        for _, idx_2, name_2 in list_len_idx_name_2:
            if idx_2 == -1:
                array_outlier_2 = (np.ones(len(array_haplotype_2[0]))*-1).astype(int)
                array_remain_2  = array_haplotype_2
            else:
                array_outlier_2 = array_haplotype_2[idx_2]
                array_remain_2  = np.vstack([array_haplotype_2[:idx_2,:], array_haplotype_2[idx_2+1:,:]])
            rec_hap1 = np.r_[array_remain_1, [array_outlier_2]]
            rec_hap2 = np.r_[array_remain_2, [array_outlier_1]]
            if len(report_outlier(rec_hap1)) == 0 and len(report_outlier(rec_hap2)) == 0:
                print("Reconcile successfully!", name_1, name_2)
                return name_1, name_2, rec_hap1, rec_hap2
    print("Reconcile Fail!!!!!")
    return None, None, None, None
    

def read_fasta(fn_fasta):
    dict_fasta = {}
    name = ""
    seq =  ""
    with open(fn_fasta) as f:
        for line in f:
            if line[0] == '>':
                if name != "":
                    dict_fasta[name] = seq
                name = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
        dict_fasta[name] = seq
    return dict_fasta





def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv1', '--gene_csv_H1', required=True, help='the IGH gene csv report for the first haplotype')
    parser.add_argument('-csv2', '--gene_csv_H2', required=True, help='the IGH gene csv report for the second haplotype')
    parser.add_argument('-f1', '--fasta_H1', help='the H1 fasta file for reconcile')
    parser.add_argument('-f2', '--fasta_H2', help='the H2 fasta file for reconcile')
    args = parser.parse_args(arguments)
    
    csv_H1 = args.gene_csv_H1
    csv_H2 = args.gene_csv_H2
    fa_H1  = args.fasta_H1
    fa_H2  = args.fasta_H2

    header, dict_contigs_genes_1, = parse_csv(csv_H1)
    header, dict_contigs_genes_2, = parse_csv(csv_H2)

    deletion_range, complex_range = sv_range(header)
    
    #print("\t".join(["","del_1", "del_5", "del_6", "del_7", "del_8", "complex"]))

    list_csv1_haplotype = []
    for contig_name, list_genes in sorted(dict_contigs_genes_1.items()):
        haplotype = report_haplotype(list_genes, deletion_range, complex_range)
        #print(contig_name, "\t".join([str(ele) for ele in haplotype]))
        list_csv1_haplotype.append(haplotype)

    list_csv2_haplotype = []
    for contig_name, list_genes in sorted(dict_contigs_genes_2.items()):
        haplotype = report_haplotype(list_genes, deletion_range, complex_range)
        #print(contig_name, "\t".join([str(ele) for ele in haplotype]))
        list_csv2_haplotype.append(haplotype)

    flag_1, list_consensus_1 = check_consistent(list_csv1_haplotype)
    flag_2, list_consensus_2 = check_consistent(list_csv2_haplotype)
    print("hap1", flag_1)
    print("hap2", flag_2)

    # if provide fasta, do reconcile
    if not check_consistent(list_csv1_haplotype)[0] or not check_consistent(list_csv2_haplotype)[0]:
        name_1, name_2, rec_hap1, rec_hap2 = reconcile(list_csv1_haplotype, list_csv2_haplotype, dict_contigs_genes_1, dict_contigs_genes_2)
        if (name_1 or name_2) and fa_H1 and fa_H2:
            dict_fa_H1 = read_fasta(fa_H1)
            dict_fa_H2 = read_fasta(fa_H2)
            f1 = open(fa_H1 + ".rec.fa", "w")
            for name, seq in dict_fa_H1.items():
                if name == name_1:
                    continue
                else:
                    f1.write('>'+name+'\n')
                    f1.write(seq+'\n')
            if name_2 != -1:
                f1.write('>'+name_2+'\n')
                f1.write(dict_fa_H2[name_2])
            f1.close()
            f2 = open(fa_H2 + ".rec.fa", "w")
            for name, seq in dict_fa_H2.items():
                if name == name_2:
                    continue
                else:
                    f2.write('>'+name+'\n')
                    f2.write(seq+'\n')
            if name_1 != -1:
                f2.write('>'+name_1+'\n')
                f2.write(dict_fa_H1[name_1])
            f2.close()
            
            # report the reconcile haplotypes
            f1 = open(fa_H1 + ".rec.csv", "w")
            f1.write(','.join([str(ele) for ele in check_consistent(rec_hap1)[1] ]))
            f1.close()
            f2 = open(fa_H2 + ".rec.csv", "w")
            f2.write(','.join([str(ele) for ele in check_consistent(rec_hap2)[1] ]))
            f2.close()
    else: # no conflict, nothing to reconcile
        if fa_H1 and fa_H2:
            subprocess.run(' '.join(["cp", fa_H1, fa_H1+".rec.fa"]), shell=True)
            subprocess.run(' '.join(["cp", fa_H2, fa_H2+".rec.fa"]), shell=True)
            f1 = open(fa_H1 + ".rec.csv", "w")
            f1.write(','.join([str(ele) for ele in list_consensus_1 ]))
            f1.close()
            f2 = open(fa_H2 + ".rec.csv", "w")
            f2.write(','.join([str(ele) for ele in list_consensus_2 ]))
            f2.close()




if __name__ == "__main__":
    main()
