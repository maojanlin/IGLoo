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


def output_bed(fo_bed, dict_read):
    if fo_bed:
        fo = open(fo_bed, 'w')
    
    for read_name, seq in dict_read.items():
        start_pos = -1
        flag_upper = False
        for idx, base in enumerate(seq):
            if base.islower() or base == 'N':
                if flag_upper:
                    if start_pos != -1:
                        write_line = read_name + '\t' + str(start_pos) + '\t' + str(idx) + '\t' + str(idx - start_pos)
                        if fo_bed:
                            fo.write(write_line + '\n')
                        else:
                            print(write_line)
                    flag_upper = False
            else:
                if flag_upper == False:
                    flag_upper = True
                    start_pos = idx
        if base.isupper():
            write_line = read_name + '\t' + str(start_pos) + '\t' + str(idx) + '\t' + str(idx - start_pos)
            if fo_bed:
                fo.write(write_line + '\n')
            else:
                print(write_line)
    if fo_bed:
        fo.close()




def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', '--input_fasta', help='targeted fasta with upper case and lower case.')
    parser.add_argument('-out',   '--output_bed',  help='the output bed file showing the upper case regions.')
    args = parser.parse_args(arguments)

    dict_read = read_fasta(args.input_fasta)
    fo_bed = args.output_bed
    output_bed(fo_bed, dict_read)



if __name__ == "__main__":
    main()

