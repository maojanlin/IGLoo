import argparse


def filter_unrecombine(input_text):
    flag_recomb = True
    for line in input_text:
        if line[:2] == "IG":
            continue
        elif line[:12] == "Unrecombined":
            flag_recomb = False
        elif line.isnumeric():
            continue
        elif line[:5] == "=====":
            flag_recomb = not flag_recomb
        else:
            if flag_recomb:
                print(line)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', '--input_text', nargs='+', default=[], help='the output of analyze_recombine.py with -name option')
    args = parser.parse_args()
    
    input_text = args.input_text
    filter_unrecombine(input_text)
    
