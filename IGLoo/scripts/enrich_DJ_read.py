import argparse
from scripts.analyze_pacbio_refs import read_fasta




def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', '--fasta', required=True, help='the fasta file for enrichment')
    parser.add_argument('-out',   '--output_fasta',  help='the output report path')
    parser.add_argument('-t',     '--enrich_time', type=int, default=1, help='number of enrichment' )
    args = parser.parse_args(arguments)

    fn_fasta = args.fasta
    fn_out   = args.output_fasta
    if fn_out == None:
        fn_out = fn_fasta + '.enrich_DJ.fa'
    enrich_time = args.enrich_time + 1
    
    dict_read = read_fasta(fn_fasta)
    f = open(fn_out, 'w')
    for name, seq in dict_read.items():
        if ("D_read" in name) or ("J_read" in name) or ("Unrecombined" in name):
            for idx in range(enrich_time):
                f.write('>' + name + '/enrich' + str(idx) + '\n')
                f.write(seq + '\n')
        else:
            f.write('>' + name + '\n')
            f.write(seq + '\n')
    f.close()




if __name__ == "__main__":
    main()

