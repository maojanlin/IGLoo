import argparse
    


def get_sample_bad(fn_bad, target_id):
    dict_sample_read = {}
    f = open(fn_bad, 'r')
    for line in f:
        fields = line.strip().split(',')
        sample_id = fields[0]
        if dict_sample_read.get(sample_id):
            dict_sample_read[sample_id].add(fields[1])
        else:
            dict_sample_read[sample_id] = {fields[1]}

    if dict_sample_read.get(target_id):
        for seq in dict_sample_read[target_id]:
            print(seq)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-bl', '--bad_log', help='the log record all the bad split mapping reads.')
    parser.add_argument('-sl', '--sample',  help='sample id we are interested in.')
    args = parser.parse_args()
    
    fn_bad = args.bad_log
    sample = args.sample
    get_sample_bad(fn_bad, sample)
    

