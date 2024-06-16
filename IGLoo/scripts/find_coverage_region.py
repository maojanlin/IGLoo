import argparse



def read_depth_log(rd_rpt):
    list_pair = []
    with open(rd_rpt) as f:
        prev_name = ""
        prev_pos = '0'
        prev_chr = ""
        flag_0 = True
        for line in f:
            chr_name, pos, depth = line.split()
            depth = int(depth)
            if chr_name != prev_name:
                if flag_0 == False:
                    list_pair.append((prev_chr,bg,prev_pos, str(int(prev_pos)-int(bg))))
                flag_0 = True
                prev_name = chr_name

            if depth > 0 and flag_0:
                bg = pos
                flag_0 = False
            elif depth == 0 and flag_0 == False:
                list_pair.append((chr_name,bg,pos, str(int(pos)-int(bg))))
                flag_0 = True
            prev_pos = pos
            prev_chr = chr_name
    if flag_0 == False:
        list_pair.append((chr_name,bg,pos, str(int(pos)-int(bg))))
    return list_pair


def output_bed(fo_bed, list_region):
    if fo_bed:
        fo = open(fo_bed, 'w')
    for ele in list_region:
        if fo_bed:
            fo.write('\t'.join(ele) + '\n')
        else:
            print('\t'.join(ele))
    if fo_bed:
        fo.close()



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-rd',  '--read_depth_log', help='the read depth log of the realignment to the region.')
    parser.add_argument('-out', '--output_bed',  help='the output bed file showing the continuous coverage region.')
    args = parser.parse_args(arguments)

    rd_rpt = args.read_depth_log
    list_region = read_depth_log(rd_rpt)
    fo_bed = args.output_bed
    output_bed(fo_bed, list_region)



if __name__ == "__main__":
    main()
