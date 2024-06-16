samtools fastq extract_hprc_bam/NA18906.hprc.long.IGH.bam > sequence_files/NA18906.hprc.long.IGH.fastq
grep -w -A 3 -f NA18906.bad.log sequence_files/NA18906.hprc.long.IGH.fastq --no-group-separator > sequence_files/NA18906.hprc.bad.fastq

hg19="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa"
minimap2 -ax map-hifi ${hg19} sequence_files/NA18906.hprc.bad.fastq > hg37_bam/NA18906.hprc.bad.bam
samtools sort hg37_bam/NA18906.hprc.bad.bam > hg37_bam/NA18906.hprc.bad.sorted.bam
samtools index hg37_bam/NA18906.hprc.bad.sorted.bam

python3 analyze_recombine.py -bed materials/gene_annotations/hg19_IGH.bed -bam hg37_bam/NA18906.hprc.bad.sorted.bam
