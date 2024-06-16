dir="example/ReAsm_out/"
sample_ID="HG005"

python3 IGLoo/scripts/ig_SV_typing.py \
    -csv1 ${dir}/asm_annotate/${sample_ID}.contig_gene.1.csv \
    -csv2 ${dir}/asm_annotate/${sample_ID}.contig_gene.2.csv \
    -f1 ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap1.fa \
    -f2 ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap2.fa

ref_base='IGLoo/materials/personalized_ref/referece_base.fa'
ref_del='IGLoo/materials/personalized_ref/referece_del.fa'

python3 IGLoo/scripts/merge_personal_ref_2hap.py \
    -base ${ref_base} \
    -alt  ${ref_del}  \
    -csv1  ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap1.fa.rec.csv \
    -csv2  ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap2.fa.rec.csv \
    -out1  ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap1.ref.fa \
    -out2  ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap2.ref.fa


echo "========  MaSuRCA ==============================="
mkdir -p ${dir}/MaSuRCA_${sample_ID}/
cd ${dir}/MaSuRCA_${sample_ID}/
cp ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap1.ref.fa ./
cp ${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap2.ref.fa ./

Ref_pat="${sample_ID}.IGH.asm.hap1.ref.fa"
Ref_mat="${sample_ID}.IGH.asm.hap2.ref.fa"
contigs_1="${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap1.fa"
contigs_2="${dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.hap2.fa"

# Arrange the contigs
chromosome_scaffolder.sh -r ${Ref_pat} -c 50000 -i 99 -m 250000 -q ${contigs_1} -t 8 -nb
chromosome_scaffolder.sh -r ${Ref_mat} -c 50000 -i 99 -m 250000 -q ${contigs_2} -t 8 -nb
arranged_contigs_1=${Ref_pat}.$(basename ${contigs_1}).split.reconciled.fa
arranged_contigs_2=${Ref_mat}.$(basename ${contigs_2}).split.reconciled.fa

# split alignments
splitScaffoldsAtNs.sh ${arranged_contigs_1} 1 > ${arranged_contigs_1}.split
splitScaffoldsAtNs.sh ${arranged_contigs_2} 1 > ${arranged_contigs_2}.split

# draft polish H1
mkdir -p draft_polish_H1 && cd draft_polish_H1
final_polish.sh 14 ../${Ref_pat} ../${Ref_pat} ../${arranged_contigs_1}
cd -

mkdir -p draft_polish_H2 && cd draft_polish_H2
final_polish.sh 14 ../${Ref_mat} ../${Ref_mat} ../${arranged_contigs_2}
cd -

python3 ../merge_haplotypes.py -f1 draft_polish_H1/14.dir/14.all.polished.fa -f2 draft_polish_H2/14.dir/14.all.polished.fa -out draft_polish.fa
minimap2 -ax map-pb draft_polish.fa ../../processed_fasta/${sample}.split.enrich.fa | samtools sort -o ${sample}.draft.realign.bam
samtools index ${sample}.draft.realign.bam


python3 ../separate_reads.py -bam ${sample}.draft.realign.bam -fasta ../../processed_fasta/${sample}.split.enrich.fa -out ${sample}.separate.read

export PYTHONPATH=/home/mlin77/lib/python3.9
cd draft_polish_H1
~/vast/maojan/jasper-1.0.3/bin/jasper.sh -t 16 -b 800000000 -a 14.dir/14.all.polished.fa -r ../${sample}.separate.read.H1.fa -k 25 -p 3
#~/vast/maojan/jasper-1.0.3/bin/jasper.sh -t 16 -b 800000000 -a 14.dir/14.all.polished.fa -r ../../../processed_fasta/${sample}.split.enrich.fa -k 25 -p 2
cd -

cd draft_polish_H2
~/vast/maojan/jasper-1.0.3/bin/jasper.sh -t 16 -b 800000000 -a 14.dir/14.all.polished.fa -r ../${sample}.separate.read.H2.fa -k 25 -p 3
#~/vast/maojan/jasper-1.0.3/bin/jasper.sh -t 16 -b 800000000 -a 14.dir/14.all.polished.fa -r ../../../processed_fasta/${sample}.split.enrich.fa -k 25 -p 2
cd -

#python3 ../merge_haplotypes.py -f1 draft_polish_H1/14.dir/14.all.polished.fa -f2 draft_polish_H2/14.dir/14.all.polished.fa -out draft_polish.fa
python3 ../merge_haplotypes.py -f1 draft_polish_H1/14.all.polished.fa.polished.fasta -f2 draft_polish_H2/14.all.polished.fa.polished.fasta -out final_polish.fa
minimap2 -ax map-pb final_polish.fa ../../processed_fasta/${sample}.split.enrich.fa | samtools sort -o ${sample}.final.realign.bam
samtools index ${sample}.final.realign.bam

python3 ../find_upperCase.py -fasta final_polish.fa -out final_polish.upperCase.bed 

samtools depth ${sample}.final.realign.bam -a -g 0x100 -J > ${sample}.realign.rd.log
python3 ../find_coverage_region.py -rd ${sample}.realign.rd.log -out ${sample}.realign.rd.bed
python3 ../crop_reassembly.py -rd ${sample}.realign.rd.bed -up final_polish.upperCase.bed -fa final_polish.fa -out final_polish.crop
cp final_polish.crop.1.fa ../crop_assembly/${sample}.crop.1.fa
cp final_polish.crop.2.fa ../crop_assembly/${sample}.crop.2.fa
cd ..
python3 ~/scr4_blangme2/maojan/IGLoo/IGLoo_asm.py -rd IGLoo_dir -id ${sample} -a1 ./crop_assembly/${sample}.crop.1.fa -a2 ./crop_assembly/${sample}.crop.2.fa

