import argparse

def parse_yak(yak_file):
    dict_yak = {}
    with open(yak_file, "r") as f:
        for line in f:
            fields = line.strip().split()
            dict_yak[fields[0]] = fields[1]
    return dict_yak

def parse_fasta(fasta_file, dict_yak, output):
    f1 = open(output + '.H1.fasta', 'w')
    f2 = open(output + '.H2.fasta', 'w')
    with open(fasta_file, "r") as f:
        flag_f1 = False
        flag_f2 = False
        for line in f:
            if line.startswith(">"):
                header = line.strip().split()[0][1:]
                if dict_yak[header] == 'p':
                    flag_f1 = True
                    flag_f2 = False
                    f1.write(line)
                elif dict_yak[header] == 'm':
                    flag_f1 = False
                    flag_f2 = True
                    f2.write(line)
                else:
                    flag_f1 = True
                    flag_f2 = True
                    f1.write(line)
                    f2.write(line)
            else:
                if flag_f1:
                    f1.write(line)
                if flag_f2:
                    f2.write(line)
    f1.close()
    f2.close()

def main(arguments=None):
    parser = argparse.ArgumentParser(description="Split read from yak")
    parser.add_argument("--input_fasta", type=str, help="Input fasta file")
    parser.add_argument("--input_yak", type=str, help="Input yak file")
    parser.add_argument("--output", type=str, help="Output header of the H1 and H2 read files")
    args = parser.parse_args(arguments)

    dict_yak = parse_yak(args.input_yak)
    parse_fasta(args.input_fasta, dict_yak, args.output)


if __name__ == "__main__":
    main()

