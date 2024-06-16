# parameters
out_dir="./example/read_out/"

# more parameters
sample_ID="HG005"
sample_fasta="example/read_out/HG005.hprc.IGH.fastq"
sample_bam_0="example/HG005.hprc.IGH.bam"
sample_bam_1=${out_dir}${sample_ID}".grch38.bam"
sample_bam_2=${out_dir}${sample_ID}".grch37.bam"
sample_bam_3=${out_dir}${sample_ID}".chm13.bam"

ref_1="/home/mlin77/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ref_2="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa" 
ref_3="/home/mlin77/data_blangme2/fasta/chm13_v2.0/chm13v2.0.fasta"

bed_1="IGLoo/materials/gene_annotations/GRCh38/grch38_IGH.bed"
bed_2="IGLoo/materials/gene_annotations/hg19_IGH.bed"
bed_3="IGLoo/materials/gene_annotations/chm13_IGH.bed"

## Do alignment
minimap2 -ax map-hifi ${ref_1} "example/read_out/HG005.IGH.fastq" | samtools sort > ${sample_bam_1}
samtools index ${sample_bam_1}
minimap2 -ax map-hifi ${ref_2} "example/read_out/HG005.IGH.fastq" | samtools sort > ${sample_bam_2}
samtools index ${sample_bam_2}
minimap2 -ax map-hifi ${ref_3} "example/read_out/HG005.IGH.fastq" | samtools sort > ${sample_bam_3}
samtools index ${sample_bam_3}

samtools fasta ${sample_bam_0} > ${sample_fasta}


mkdir ${out_dir}processed_fasta
mkdir ${out_dir}pc_report

echo "[Deplete ReAssembly] Processing ${sample_ID}..."
python3 IGLoo/scripts/analyze_pacbio_refs.py -lbed ${bed_1} \
                                     ${bed_1} \
                                     ${bed_2} \
                                     ${bed_3} \
                               -lbam ${sample_bam_0} \
                                     ${sample_bam_1} \
                                     ${sample_bam_2} \
                                     ${sample_bam_3} \
                               -fasta ${sample_fasta} \
                               -out_fa ${out_dir}/processed_fasta/${sample_ID}.split.fa \
                               -out_dt ${out_dir}/pc_report/${sample_ID}.split.detail.rpt \
                               -out    ${out_dir}/pc_report/${sample_ID}.split.rpt

python3 IGLoo/scripts/enrich_DJ_read.py -fasta ${out_dir}/processed_fasta/${sample_ID}.split.fa \
                                        -out   ${out_dir}/processed_fasta/${sample_ID}.split.enrich.fa


#parent_1_yak="parent_data_clipped/${sample_ID}/pat.yak"
#parent_2_yak="parent_data_clipped/${sample_ID}/mat.yak"

#"""
#python3 IGLoo/scripts/analyze_pacbio_refs.py 
#
#-h
#usage: analyze_pacbio_refs.py [-h] -lbed LIST_ANNOTATION_BED [LIST_ANNOTATION_BED ...] -lbam LIST_ALIGN_BAM [LIST_ALIGN_BAM ...] -fasta FASTA
#                              [-out OUTPUT_REPORT] [-out_dt OUTPUT_DETAIL] -out_fa OUTPUT_FASTA
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -lbed LIST_ANNOTATION_BED [LIST_ANNOTATION_BED ...], --list_annotation_bed LIST_ANNOTATION_BED [LIST_ANNOTATION_BED ...]
#                        the IGH annotations of the reference genomes.
#  -lbam LIST_ALIGN_BAM [LIST_ALIGN_BAM ...], --list_align_bam LIST_ALIGN_BAM [LIST_ALIGN_BAM ...]
#                        the alignment bam files of the IG reads to the reference genome.
#  -fasta FASTA, --fasta FASTA
#                        the fasta file of the bams
#  -out OUTPUT_REPORT, --output_report OUTPUT_REPORT
#                        the path of output report
#  -out_dt OUTPUT_DETAIL, --output_detail OUTPUT_DETAIL
#                        the path of detailed report
#  -out_fa OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
#                        the output split read fasta
#(base) [mlin77@devlangmead1 17:32:03 IGLoo]$ python3 IGLoo/IGLoo.py -s HG005 -o example/read_out/ -b example/HG005.hprc.IGH.bam -lb IGLoo/materials/gene_annotations/GRCh38/grch38_IGH.bed IGLoo/materials/gene_annotations/hg19_IGH.bed IGLoo/materials/gene_annotations/chm13_IGH.bed -lr /home/mlin77/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /home/mlin77/data_blangme2/fasta/grch37/hg19.fa /home/mlin77/data_blangme2/fasta/chm13_v2.0/chm13v2.0.fasta
#"""
