#wgs_bam="/home/mlin77/data_blangme2/reads/HPRC/HG00733/HG00733_aligned_GRCh38_winnowmap.sorted.bam"
#batch_id="HG00733.hprc"
#wgs_bam="/home/mlin77/data_blangme2/reads/HPRC/NA19240/NA19240_aligned_GRCh38_winnowmap.sorted.bam"
#batch_id="NA19240.hprc"
wgs_bam="/home/mlin77/data_blangme2/reads/primary_cell/SRS1156159.grch38.bam"
batch_id="SAMN04251426.primary"
out_dir="extract_bam"
IGH_bed="materials/IGH_grch38_extend.bed"
IGL_bed="materials/IGL_grch38_extend.bed"
IGK_bed="materials/IGK_grch38_extend.bed"

samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
samtools index ${out_dir}/${batch_id}.IGH.bam
samtools index ${out_dir}/${batch_id}.IGK.bam
samtools index ${out_dir}/${batch_id}.IGL.bam

#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 
