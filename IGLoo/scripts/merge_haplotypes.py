import argparse


def read_fasta(fn_fasta):
    dict_read = {}
    name = ""
    seq  = ""
    with open(fn_fasta) as f:
        for line in f:
            if line[0] == '>':
                if name != "":
                    dict_read[name] = seq
                name = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
    dict_read[name] = seq
    return dict_read


def output_merge(fo_merge, dict_read_1, dict_read_2):
    fo = open(fo_merge, 'w')
    for name, seq in dict_read_1.items():
        fo.write('>' + name + '_H1\n')
        fo.write(seq + '\n')
    for name, seq in dict_read_2.items():
        fo.write('>' + name + '_H2\n')
        fo.write(seq + '\n')
    fo.close()



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', '--input_fasta_H1', help='draft H1 fasta.')
    parser.add_argument('-f2', '--input_fasta_H2', help='draft H2 fasta.')
    parser.add_argument('-out', '--output_merge_fasta',  help='the output merged fasta file.')
    args = parser.parse_args(arguments)

    dict_read_1 = read_fasta(args.input_fasta_H1)
    dict_read_2 = read_fasta(args.input_fasta_H2)
    output_merge(args.output_merge_fasta, dict_read_1, dict_read_2)




if __name__ == "__main__":
    main()
