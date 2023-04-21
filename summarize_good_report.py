import argparse


def summarize_report(fn_report, dict_full, dict_VDJ):
    f = open(fn_report)
    for line in f:
        combination, occurence = line.strip().split('\t')
        combination = tuple(combination.split('---'))
        comb_VDJ = []
        for gene in combination:
            if gene[3] in {'V', 'D', 'J', 'e'}:
                comb_VDJ.append(gene)
        comb_VDJ = tuple(comb_VDJ)
        occurence = int(occurence)

        if dict_full.get(combination):
            dict_full[combination] += occurence
        else:
            dict_full[combination] = occurence
        if dict_VDJ.get(comb_VDJ):
            dict_VDJ[comb_VDJ] += occurence
        else:
            dict_VDJ[comb_VDJ] = occurence
    f.close()






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-lr', '--list_report', nargs='+', required=True, help='the list of the good reports')
    parser.add_argument('-bad', '--bad_file', help="bad file report.")
    parser.add_argument('-out', '--output_file', help="specify for the final good report result.")
    parser.add_argument('-const', '--constant_gene', action='store_true', help='consider constant gene in classification')
    args = parser.parse_args()
    
    list_report = args.list_report
    bad_file    = args.bad_file
    output_file = args.output_file
    flag_const  = args.constant_gene

    dict_full = {}
    dict_VDJ = {}
    for fn_report in list_report:
        summarize_report(fn_report, dict_full, dict_VDJ)
    
    unknown_num = 0
    if bad_file:
        f = open(bad_file)
        for line in f:
            unknown_num += 1
        f.close()
    if unknown_num > 0:
        dict_full[("Unknown",)] = unknown_num
        dict_VDJ [("Unknown",)] = unknown_num
    
    if output_file:
        f = open(output_file, 'w')
        if flag_const:
            for combination, number in sorted(dict_full.items()):
                f.write('---'.join(combination) + '\t' + str(number) + '\n')
        else:
            for combination, number in sorted(dict_VDJ.items()):
                f.write('---'.join(combination) + '\t' + str(number) + '\n')
        f.close()
    total_num = sum(dict_VDJ.values())
    max_num   = max(dict_VDJ.values())
    if dict_VDJ.get(("Unrecombined",)):
        unrecombine_num = dict_VDJ[("Unrecombined",)]
    else:
        unrecombine_num = 0
    print(len(dict_full), len(dict_VDJ), total_num, max_num, unrecombine_num, sep='\t')
