#wgs_bam="/home/mlin77/data_blangme2/reads/HPRC/HG00733/HG00733_aligned_GRCh38_winnowmap.sorted.bam"
#batch_id="HG00733.hprc"
wgs_bam="/home/mlin77/data_blangme2/reads/HPRC/NA18906_aligned_GRCh38_winnowmap.sorted.bam"
batch_id="NA18906.hprc.long"
out_dir="extract_hprc_bam"
IGH_bed="materials/IGH_grch38_extend.bed"
#IGL_bed="materials/IGL_grch38_extend.bed"
##IGK_bed="materials/IGK_grch38_extend.bed"

#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
samtools view -h ${wgs_bam} "chr14:105486937-106979845" -@ 8 -o ${out_dir}/${batch_id}.IGH.bam
##samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
##samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
samtools index ${out_dir}/${batch_id}.IGH.bam
##samtools index ${out_dir}/${batch_id}.IGK.bam
##samtools index ${out_dir}/${batch_id}.IGL.bam

#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 
