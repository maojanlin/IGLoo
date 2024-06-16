for sample_id in HG002 \
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
do
    echo ${sample_id}
    IGHC="/scratch4/blangme2/maojan/AIRRssembly/constant_bam/IGHC.fasta"
    dir_asm="/home/mlin77/data_blangme2/fasta/hprc-yr1"
    
    minimap2 -ax map-hifi ${dir_asm}/${sample_id}.1.fa ${IGHC} > constant_bam/hprc.${sample_id}.1.IGHC.bam
    samtools sort constant_bam/hprc.${sample_id}.1.IGHC.bam > constant_bam/hprc.${sample_id}.1.IGHC.sorted.bam
    samtools index constant_bam/hprc.${sample_id}.1.IGHC.sorted.bam
    minimap2 -ax map-hifi ${dir_asm}/${sample_id}.2.fa ${IGHC} > constant_bam/hprc.${sample_id}.2.IGHC.bam
    samtools sort constant_bam/hprc.${sample_id}.2.IGHC.bam > constant_bam/hprc.${sample_id}.2.IGHC.sorted.bam
    samtools index constant_bam/hprc.${sample_id}.2.IGHC.sorted.bam
    rm constant_bam/hprc.${sample_id}.1.IGHC.bam
    rm constant_bam/hprc.${sample_id}.2.IGHC.bam
done

