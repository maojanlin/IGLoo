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


def read_bed(fn_bed):
    dict_region = {}
    f = open(fn_bed)
    for line in f:
        fields = line.strip().split()
        if dict_region.get(fields[0]):
            dict_region[fields[0]].append(fields[1:])
        else:
            dict_region[fields[0]] = [fields[1:]]
    f.close()
    return dict_region


def find_overlap(dict_rd, dict_up):
    dict_crop_region = {}
    for contig in sorted(dict_rd.keys()):
        dict_crop_region[contig] = []
        list_up_region = []
        for bg, ed, length in dict_up[contig]:
            if int(length) > 100:
                list_up_region.append((int(bg), int(ed)))
        list_rd_region = [(int(ele[0]), int(ele[1])) for ele in dict_rd[contig]]
        for rd_bg, rd_ed in list_rd_region:
            for up_bg, up_ed in list_up_region:
                if (rd_bg <= up_bg and rd_ed >= up_bg) or (rd_bg <= up_ed and rd_ed >= up_ed): # include the upperCase region
                    dict_crop_region[contig].append((rd_bg, rd_ed))
                    break
    print(dict_crop_region)
    return dict_crop_region


def output_fasta(dict_crop_region, dict_fasta, fo_fasta):
    for idx, contig in enumerate(["chr14_H1", "chr14_H2"]):
        fo = open(fo_fasta+"."+str(idx+1)+".fa", "w")
        for idy, region in enumerate(dict_crop_region[contig]):
            bg, ed = region
            fo.write(">IGH_contig_"+str(idy+1)+"\n")
            fo.write(dict_fasta[contig][bg:ed] + "\n")
        fo.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-rd',  '--read_depth_bed', required=True, help='the read depth log of the realignment to the region.')
    parser.add_argument('-up',  '--upperCase_bed', required=True, help='the upper case bed report from the reassembly')
    parser.add_argument('-fa',  '--reassembly_fa', required=True, help='the reassembled fasta')
    parser.add_argument('-out', '--output_fa',  help='the output fasta of the two haplotypes.')
    args = parser.parse_args()

    rd_bed = args.read_depth_bed
    up_bed = args.upperCase_bed
    fn_fasta = args.reassembly_fa
    fo_fasta = args.output_fa

    dict_fasta = read_fasta(fn_fasta)
    dict_rd = read_bed(rd_bed)
    dict_up = read_bed(up_bed)

    dict_crop_region = find_overlap(dict_rd, dict_up)
    output_fasta(dict_crop_region, dict_fasta, fo_fasta)

