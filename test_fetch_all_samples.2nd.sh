###### Get the fastq file ######
grep -w -A 3 -f candidate_2_alt_ref_2nd.log sequence_files/all_sample.bad.fastq --no-group-separator > sequence_files/all_sample.bad.2nd.fastq
###### Get the fastq file ######

chm13="/home/mlin77/data_blangme2/fasta/chm13_v2.0/chm13v2.0.fasta"
minimap2 -ax map-hifi ${chm13} sequence_files/all_sample.bad.2nd.fastq > chm13_bam/all_sample.hprc.bad.2nd.bam
samtools sort chm13_bam/all_sample.hprc.bad.2nd.bam > chm13_bam/all_sample.hprc.bad.2nd.sorted.bam
samtools index chm13_bam/all_sample.hprc.bad.2nd.sorted.bam

python3 analyze_recombine.py -bed materials/gene_annotations/chm13_IGH.bed -bam chm13_bam/all_sample.hprc.bad.2nd.sorted.bam
