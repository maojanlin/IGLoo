crop_hprc(){
    batch_id=$1
    wgs_bam="/home/mlin77/data_blangme2/reads/HPRC/${batch_id}_aligned_GRCh38_winnowmap.sorted.bam"
    out_dir="extract_hprc_bam"
    IGH_bed="materials/IGH_grch38_extend.bed"
    
    samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.hprc.IGH.bam
    samtools index ${out_dir}/${batch_id}.hprc.IGH.bam
}
export -f crop_hprc
parallel crop_hprc ::: \
    HG002   \
    HG00438 \
    HG005   \
    HG00621 \
    HG00673 \
    HG00733 \
    HG00735 \
    HG00741 \
    HG01071 \
    HG01106 \
    HG01109 \
    HG01123 \
    HG01175 \
    HG01243 \
    HG01258 \
    HG01358 \
    HG01361 \
    HG01891 \
    HG01928 \
    HG01952 \
    HG01978 \
    HG02055 \
    HG02080 \
    HG02109 \
    HG02145 \
    HG02148 \
    HG02257 \
    HG02486 \
    HG02559 \
    HG02572 \
    HG02622 \
    HG02630 \
    HG02717 \
    HG02723 \
    HG02818 \
    HG02886 \
    HG03098 \
    HG03453 \
    HG03486 \
    HG03492 \
    HG03516 \
    HG03540 \
    HG03579 \
    NA18906 \
    NA19240 \
    NA20129 \
    NA21309
