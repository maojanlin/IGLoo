# The goal of this script is to map IGH reads to all the three reference genomes for aggregate analysis
map_to_references() {
    # parameters
    out_dir="IGH_mappings/"
    mkdir -p ${out_dir}
    
    reference_1="/home/mlin77/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    reference_2="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa"
    reference_3="/home/mlin77/data_blangme2/fasta/chm13_v2.0/chm13v2.0.fasta"

    sample_ID=${1}
    sample_fasta="depleted_fasta/${sample_ID}.hprc.IGH.fasta"

    echo "Map $sample_ID to GRCh38"
    minimap2 -ax map-hifi ${reference_1} ${sample_fasta} | samtools sort -o ${out_dir}/${sample_ID}.IGH.grch38.bam
    samtools index ${out_dir}/${sample_ID}.IGH.grch38.bam
    echo "Map $sample_ID to GRCh37"
    minimap2 -ax map-hifi ${reference_2} ${sample_fasta} | samtools sort -o ${out_dir}/${sample_ID}.IGH.grch37.bam
    samtools index ${out_dir}/${sample_ID}.IGH.grch37.bam
    echo "Map $sample_ID to CHM13"
    minimap2 -ax map-hifi ${reference_3} ${sample_fasta} | samtools sort -o ${out_dir}/${sample_ID}.IGH.chm13.bam
    samtools index ${out_dir}/${sample_ID}.IGH.chm13.bam
}
export -f map_to_references 
parallel map_to_references ::: \
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

