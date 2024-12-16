#Force polish a genome with samtools mpileup

import argparse
from collections import defaultdict

def parse_bases(bases: str) -> dict:
    """Parse mpileup base string and return counts of each base type.
    
    Args:
        bases: String from mpileup containing base information
              '.' or ',' = match to reference
              'ACGT' or 'acgt' = mismatch base
              '+Nx' = insertion of N bases
              '-Nx' = deletion of N bases (preceded by . or ,)
              
    Returns:
        Dictionary with counts of each base type, including deletions
    """
    parsed_bases = defaultdict(int)
    i = 0
    while i < len(bases):
        # Skip start/end of read markers
        if bases[i] == '^' and i + 1 < len(bases):
            i += 2  # Skip quality score after ^
            continue
        if bases[i] == '$':
            i += 1
            continue
            
        # Handle deletions (must check for preceding . or ,)
        if i + 1 < len(bases) and bases[i] in '.,' and bases[i+1] == '-':
            i += 2  # Skip the . and -
            # Get length of deletion
            del_len = 0
            while i < len(bases) and bases[i].isdigit():
                del_len = del_len * 10 + int(bases[i])
                i += 1
            # Get the deleted sequence
            del_seq = bases[i:i+del_len]
            parsed_bases[f'-{del_seq.upper()}'] += 1
            i += del_len
            continue
            
        # Handle insertions
        if bases[i] == '+':
            i += 1
            # Get length of insertion
            ins_len = 0
            while i < len(bases) and bases[i].isdigit():
                ins_len = ins_len * 10 + int(bases[i])
                i += 1
            # Get the inserted sequence
            ins_seq = bases[i:i+ins_len]
            parsed_bases[f'+{ins_seq.upper()}'] += 1
            i += ins_len
            continue
            
        # Count reference matches
        if bases[i] in ',.':
            parsed_bases['ref'] += 1
            
        # Count alternative bases
        elif bases[i].upper() in 'ACGT':
            parsed_bases[bases[i].upper()] += 1
            
        i += 1
        
    return parsed_bases

def read_mpileup(mpileup_file):
    dict_change = defaultdict(list)
    with open(mpileup_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom, pos, ref, depth, bases = fields[0], int(fields[1]), fields[2], int(fields[3]), fields[4]
            # if over 0.8 of the depth is one specific alt base, call it
            parsed_bases = parse_bases(bases)
            for base, count in parsed_bases.items():
                if base == 'ref':
                    continue
                if count / depth > 0.8:
                    dict_change[chrom].append((pos, base))
                    #print(f"{chrom}:{pos} {ref}->{base}")
    return dict_change

def force_polish(original_seq, list_change):
    for pos, base in list_change[::-1]:
        if original_seq[pos-1] == 'N':
            continue
        if base.startswith('-'):
            original_seq = original_seq[:pos] + original_seq[pos+len(base)-1:]
        elif base.startswith('+'):
            original_seq = original_seq[:pos] + base[1:] + original_seq[pos:]
        else:
            original_seq = original_seq[:pos-1] + base + original_seq[pos:]
    return original_seq

def output_polish(genome_file, dict_change, out_file):
    with open(genome_file, "r") as f, open(out_file, "w") as out:
        for line in f:
            if line.startswith(">"):
                out.write(line)
                chrom = line.strip().split(" ")[0][1:]
            else:
                polished_seq = force_polish(line.strip(), dict_change[chrom])
                out.write(polished_seq + "\n")


def main(arguments=None):
    parser = argparse.ArgumentParser(description="Force polish a genome with samtools mpileup")
    parser.add_argument("--genome", type=str, help="Path to the genome fasta file")
    parser.add_argument("--mpileup", type=str, help="Path to the mpileup file")
    parser.add_argument("--out", type=str, help="Path to the output file")
    args = parser.parse_args(arguments)

    dict_change = read_mpileup(args.mpileup)
    output_polish(args.genome, dict_change, args.out)

if __name__ == "__main__":
    main()

