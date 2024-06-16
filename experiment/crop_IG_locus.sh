#wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002.15kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam"
batch_id="HG002.15kb"
out_dir="extract_bam"
IGH_bed="materials/IGH_grch38_extend.bed"
IGL_bed="materials/IGL_grch38_extend.bed"
IGK_bed="materials/IGK_grch38_extend.bed"
#
#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
#samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
#samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
#
#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 
#
#
#wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam"
#batch_id="HG002.15kb_20kb"
#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
#samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
#samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
#
#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 
#
#
#
#wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002_GRCh38.haplotag.10x.bam"
#batch_id="HG002.11kb"
#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
#samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
#samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
#
#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 

wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002.SequelII_early.bam"
batch_id="HG002.early"
samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam

samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 

#wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002_PacBio_GRCh38.bam"
#batch_id="HG002.70x"
#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
#samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
#samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
#
#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 

#wgs_bam="/home/mlin77/data_blangme2/reads/hg002_HiFi/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam"
#batch_id="HG002.10kb"
#out_dir="extract_bam"
#IGH_bed="materials/IGH_grch37_extend.bed"
#IGL_bed="materials/IGL_grch37_extend.bed"
#IGK_bed="materials/IGK_grch37_extend.bed"
#
#samtools view -h ${wgs_bam}  -L ${IGH_bed} -o  ${out_dir}/${batch_id}.IGH.bam
#samtools view -h ${wgs_bam}  -L ${IGK_bed} -o  ${out_dir}/${batch_id}.IGK.bam
#samtools view -h ${wgs_bam}  -L ${IGL_bed} -o  ${out_dir}/${batch_id}.IGL.bam
#
#samtools fasta ${out_dir}/${batch_id}.IGH.bam > ${out_dir}/${batch_id}.IGH.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGL.bam > ${out_dir}/${batch_id}.IGL.fasta 
#samtools fasta ${out_dir}/${batch_id}.IGK.bam > ${out_dir}/${batch_id}.IGK.fasta 
