# usage: bash fetch_IGH_from_chm13.sh sample_ID input_bam output_path

IGH_locus="chr14:99739969-101161492"
sample_ID=$1
input_bam=$2
output_path=$3

mkdir -p ${output_path}
samtools view -h ${input_bam} ${IGH_locus}     -o ${output_path}/${sample_ID}.IGH_chr14.bam
samtools view -h ${input_bam} '*'              -o ${output_path}/${sample_ID}.unmapped.bam

cd ${output_path}
samtools merge ${sample_ID}.IGH.bam \
               ${sample_ID}.IGH_chr14.bam \
               ${sample_ID}.unmapped.bam
samtools fastq ${sample_ID}.IGH.bam      >  ${sample_ID}.IGH.fq

cd -
