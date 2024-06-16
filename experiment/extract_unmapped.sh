samtools view -h -f 4 ~/data_blangme2/reads/HPRC/HG03486_aligned_GRCh38_winnowmap.sorted.bam -o extract_unmapped_bam/HG03486.hprc.unmapped.bam
samtools fasta extract_unmapped_bam/HG03486.hprc.unmapped.bam > extract_unmapped_bam/HG03486.hprc.unmapped.fasta
