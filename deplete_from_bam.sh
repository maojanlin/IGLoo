sample_ID=HG03492

# making the filter list for recombined sequence
python3 analyze_recombine.py -bed materials/gene_annotations/GRCh38/grch38_IGH.bed -bam extract_hprc_bam/${sample_ID}.hprc.IGH.bam -name \
    | xargs python3 select_unrecombined.py -input > depleted_fasta/${sample_ID}.recombine.list

# filter out the recombined sequence to a fasta containing only the unrecombined sequence
samtools fasta extract_hprc_bam/${sample_ID}.hprc.IGH.bam > depleted_fasta/${sample_ID}.hprc.IGH.fasta
python3 deplete_fasta.py -fasta depleted_fasta/${sample_ID}.hprc.IGH.fasta -list depleted_fasta/${sample_ID}.recombine.list -o depleted_fasta/${sample_ID}.hprc.IGH.deleted.fasta

# hifiasm
echo "[AIRRsembly] Assemble with hifiasm..."
mkdir -p depleted_fasta/${sample_ID}_hifiasm
cd depleted_fasta/${sample_ID}_hifiasm
/home/mlin77/scr4_blangme2/maojan/hifiasm/hifiasm -o ${sample_ID}.IGH.asm -t4 -f0 ../${sample_ID}.hprc.IGH.deleted.fasta 2> ${sample_ID}.IGH.asm.log
awk '/^S/{print ">"$2;print $3}' ${sample_ID}.IGH.asm.bp.p_ctg.gfa > ${sample_ID}.IGH.asm.p_ctg.fa
cd -
