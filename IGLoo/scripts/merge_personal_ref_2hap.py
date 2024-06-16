import argparse



def read_fasta(fn_fasta):
    dict_fasta = {}
    name = ""
    seq =  ""
    with open(fn_fasta) as f:
        for line in f:
            if line[0] == '>':
                if name != "":
                    dict_fasta[name] = seq
                name = line[1:].strip().split(':')[0]
                seq = ""
            else:
                seq += line.strip()
        dict_fasta[name] = seq
    return dict_fasta


def read_csv(fn_csv):
    # Note that in input haplotype format, last bit is for the complex SV
    if fn_csv:
        with open(fn_csv) as f:
            line = f.readline()
        list_sv_hap = line.strip().split(',')
        list_sv_hap = [int(ele) for ele in list_sv_hap]
        list_sv_hap = list_sv_hap[:2] + [1 - list_sv_hap[-1]] + list_sv_hap[2:-1]
        return list_sv_hap
    else:
        print("WARNING! no csv provided, use full IGH region.")
        return [1,1,0,1,1,1,1]


def piece_together(dict_base, dict_alt, list_sv_hap):
    combined_seq = dict_base['base_seq_0']
    for idx, sv_hap in enumerate(list_sv_hap):
        if sv_hap == 0: # deletion
            combined_seq += dict_alt['sv_' + str(idx)]
        else: # hap > 1 or unknown
            combined_seq += dict_base['sv_' + str(idx)]
        combined_seq += dict_base['base_seq_' + str(idx+1)]
    return combined_seq




def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-base', '--base_fasta', required=True, help='the backbone of the IG reference genome')
    parser.add_argument('-alt',  '--alt_fasta',  required=True, help='the alternative (deletion pieces) reference genome')
    parser.add_argument('-csv1', '--SV_haplotype_1', help='the SV haplotype csv file')
    parser.add_argument('-csv2', '--SV_haplotype_2', help='the SV haplotype csv file')
    parser.add_argument('-out1', '--output_fasta_1', default="./personal_IGH.1.fa", help='the output persoanlized IG reference genome')
    parser.add_argument('-out2', '--output_fasta_2', default="./personal_IGH.2.fa", help='the output persoanlized IG reference genome')
    args = parser.parse_args(arguments)
    

    dict_base = read_fasta(args.base_fasta)
    dict_alt  = read_fasta(args.alt_fasta)
    list_sv_hap_1 = read_csv(args.SV_haplotype_1)
    list_sv_hap_2 = read_csv(args.SV_haplotype_2)

    # Try to resolve the complex region by guessing
    if list_sv_hap_1[2] == 2:
        if list_sv_hap_2[2] == 2:
            list_sv_hap_1[2] = 0
            list_sv_hap_2[2] = 1
        else:
            list_sv_hap_1[2] = 1 - list_sv_hap_2[2]
    elif list_sv_hap_2[2] == 2:
        list_sv_hap_2[2] = 1 - list_sv_hap_1[2]

    combined_seq_1 = piece_together(dict_base, dict_alt, list_sv_hap_1)
    combined_seq_2 = piece_together(dict_base, dict_alt, list_sv_hap_2)

    fo = open(args.output_fasta_1, 'w')
    fo.write('>chr14\n')
    fo.write(combined_seq_1)
    fo.close()
    fo = open(args.output_fasta_2, 'w')
    fo.write('>chr14\n')
    fo.write(combined_seq_2)
    fo.close()


if __name__ == "__main__":
    main()
