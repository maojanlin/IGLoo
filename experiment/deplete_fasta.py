import argparse
import pysam
import subprocess


def read_list(fn_list):
    set_target_read = set()
    fi = open(fn_list, 'r')
    for line in fi:
        set_target_read.add(line.strip())
    return set_target_read


def deplete_fasta(fn_fasta, fn_out, set_target_read):
    fo = open(fn_out, 'w')
    fi = open(fn_fasta, 'r')

    set_target_read.add("")
    
    name = ""
    seq = ""
    for line in fi:
        if line[0] == '>':
            if name not in set_target_read:
                fo.write('>' + name + '\n')
                fo.write(seq + '\n')
            name = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
    # last sequence
    if name not in set_target_read:
        fo.write('>' + name + '\n')
        fo.write(seq)

    fi.close()
    fo.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', '--input_fasta', help='targeted original fasta containing unwanted reads.')
    parser.add_argument('-list',  '--remove_list', help='list containing all the read name to be depleted.')
    parser.add_argument('-o', '--depleted_fasta', help='final fasta file without the unwanted reads.')
    args = parser.parse_args()


    fn_fasta = args.input_fasta
    fn_list  = args.remove_list
    fn_out   = args.depleted_fasta
    
    set_target_read = read_list(fn_list)
    deplete_fasta(fn_fasta, fn_out, set_target_read)

