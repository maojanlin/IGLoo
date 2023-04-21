import argparse


def filter_unrecombine(input_text):
    for line in input_text:
        if line[:2] == "IG":
            continue
        elif line[:12] == "Unrecombined":
            break
        elif line.isnumeric():
            continue
        else:
            print(line)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', '--input_text', nargs='+', default=[], help='the output of analyze_recombine.py with -name option')
    args = parser.parse_args()
    
    input_text = args.input_text
    filter_unrecombine(input_text)
    
