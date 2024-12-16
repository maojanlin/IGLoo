import argparse
import re


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


def mask_seq(regions, sequence):
    old_bg = 0
    for region in regions:
        bg, ed, _ = region
        bg = int(bg) -1
        sequence = sequence[:old_bg] + 'N'*(bg-old_bg) + sequence[bg:]
        #print(bg-old_bg)
        old_bg = int(ed) -1
    return sequence



def output_fasta(dict_rd, dict_fasta, fo_fasta, dict_legit):
    contig = "chr14"
    fo = open(fo_fasta, "w")
    hap_regions  = dict_rd[contig]
    hap_sequence = dict_fasta[contig]
    masked_sequence = mask_seq(hap_regions, hap_sequence)

    for idy, region in enumerate(dict_legit[contig]):
        bg, ed = region
        fo.write(">IGH_masked_"+str(idy+1)+"\n")
        #fo.write(">chr14\n")
        fo.write(masked_sequence[bg:ed] + "\n")
    fo.close()


def check_Ns(dict_fasta, dict_up):
    dict_legit = {}
    contig = "chr14"
    target_seq = dict_fasta[contig]
    target_up  = dict_up[contig]
    dict_legit[contig] = []

    segments = re.split('N+', target_seq)
    if len(segments) == 1:
        dict_legit[contig].append((0, len(target_seq)))
    else:
        for segment in segments:
            seg_start = target_seq.find(segment)
            seg_stop  = seg_start + len(segment)
            for up_info in target_up:
                up_bg, up_ed, up_len = up_info
                if int(up_len) < 1000:
                    continue
                up_bg = int(up_bg) - 1
                up_ed = int(up_ed) - 1
                if (up_bg >= seg_start and up_bg <= seg_stop) or (up_ed >= seg_start and up_ed <= seg_stop):
                    dict_legit[contig].append((seg_start, seg_stop))
                    break
        print(dict_legit)
    return dict_legit





def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-rd',  '--read_depth_bed', required=True, help='the read depth log of the realignment to the region.')
    parser.add_argument('-up',  '--upperCase_bed', required=True, help='the upper case bed report from the reassembly')
    parser.add_argument('-fa',  '--reassembly_fa', required=True, help='the reassembled fasta')
    parser.add_argument('-out', '--output_fa',  help='the output fasta of the two haplotypes.')
    args = parser.parse_args(arguments)

    rd_bed = args.read_depth_bed
    up_bed = args.upperCase_bed
    fn_fasta = args.reassembly_fa
    fo_fasta = args.output_fa

    dict_fasta = read_fasta(fn_fasta)
    dict_rd = read_bed(rd_bed)
    dict_up = read_bed(up_bed)

    dict_legit = check_Ns(dict_fasta, dict_up)
    #dict_legit = {'chr14':[(0, len(dict_fasta['chr14']))]}
    output_fasta(dict_rd, dict_fasta, fo_fasta, dict_legit)

if __name__ == "__main__":
    main()

