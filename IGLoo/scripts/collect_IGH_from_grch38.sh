# usage: bash collect_IGH_from_grch38.sh sample_ID input_bam output_path [--unmap]

IGH_locus="chr14:105486937-106979845"
IGH_alt_locus="chr14_KI270726v1_random:0-43739"
sample_ID=$1
input_bam=$2
output_path=$3
flag_unmap=$4

include_unmapped=false

# Check for --unmap option
if [[ "$flag_unmap" == "--unmap" ]]; then
    include_unmapped=true
fi

mkdir -p ${output_path}
samtools view -h ${input_bam} ${IGH_locus}     -o ${output_path}/${sample_ID}.IGH_chr14.bam
samtools view -h ${input_bam} ${IGH_alt_locus} -o ${output_path}/${sample_ID}.IGH_alt.bam
if ${include_unmapped}; then
    samtools view -h ${input_bam} '*'          -o ${output_path}/${sample_ID}.unmapped.bam
fi

cd ${output_path}
if ${include_unmapped}; then
    samtools merge ${sample_ID}.IGH.bam \
                   ${sample_ID}.IGH_chr14.bam \
                   ${sample_ID}.IGH_alt.bam \
                   ${sample_ID}.unmapped.bam
else
    samtools merge ${sample_ID}.IGH.bam \
                   ${sample_ID}.IGH_chr14.bam \
                   ${sample_ID}.IGH_alt.bam
fi
samtools fastq ${sample_ID}.IGH.bam      >  ${sample_ID}.IGH.fq

cd -
