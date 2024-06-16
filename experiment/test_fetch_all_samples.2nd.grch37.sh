hg19="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa"
minimap2 -ax map-hifi ${hg19} sequence_files/all_sample.bad.2nd.fastq > hg37_bam/all_sample.hprc.bad.2nd.bam
samtools sort hg37_bam/all_sample.hprc.bad.2nd.bam > hg37_bam/all_sample.hprc.bad.2nd.sorted.bam
samtools index hg37_bam/all_sample.hprc.bad.2nd.sorted.bam
