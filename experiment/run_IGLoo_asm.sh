[ -e asm_annotation/trio_clipping/summary.txt    ] && rm asm_annotation/trio_clipping/summary.txt
[ -e asm_annotation/trio_no_clipping/summary.txt ] && rm asm_annotation/trio_no_clipping/summary.txt
[ -e asm_annotation/clipping/summary.txt    ] && rm asm_annotation/clipping/summary.txt
[ -e asm_annotation/no_clipping/summary.txt ] && rm asm_annotation/no_clipping/summary.txt
for sample_id in \
    HG00438 \
    HG00621 \
    HG00673 \
    HG00733 \
    HG00735 \
    HG00741 \
    HG01071 \
    HG01106 \
    HG01109 \
    HG01123 \
    HG01175 \
    HG01243 \
    HG01258 \
    HG01358 \
    HG01361 \
    HG01891 \
    HG01928 \
    HG01952 \
    HG01978 \
    HG02055 \
    HG02080 \
    HG02109 \
    HG02145 \
    HG02148 \
    HG02257 \
    HG02486 \
    HG02559 \
    HG02572 \
    HG02622 \
    HG02630 \
    HG02717 \
    HG02723 \
    HG02818 \
    HG02886 \
    HG03098 \
    HG03453 \
    HG03486 \
    HG03492 \
    HG03516 \
    HG03540 \
    HG03579 \
    NA18906 \
    NA19240 \
    NA20129
do
    echo ${sample_id}
    python3 IGLoo_asm.py -rd asm_annotation/trio_clipping/ -id ${sample_id} \
        -a1 parent_data_clipped/hifiasm_result/${sample_id}.IGH.asm.hap1.fa \
        -a2 parent_data_clipped/hifiasm_result/${sample_id}.IGH.asm.hap2.fa
    python3 IGLoo_asm.py -rd asm_annotation/trio_no_clipping/ -id ${sample_id} \
        -a1 parent_data/hifiasm_result_noclip/${sample_id}.IGH.asm.hap1.fa \
        -a2 parent_data/hifiasm_result_noclip/${sample_id}.IGH.asm.hap2.fa
    python3 IGLoo_asm.py -rd asm_annotation/clipping/ -id ${sample_id} \
        -a1 clipped_fasta/hifiasm_result/${sample_id}.IGH.asm.hap1.fa \
        -a2 clipped_fasta/hifiasm_result/${sample_id}.IGH.asm.hap2.fa
    python3 IGLoo_asm.py -rd asm_annotation/no_clipping/ -id ${sample_id} \
        -a1 IGH_contig_primary/hifiasm_result/${sample_id}.IGH.p_ctg.fa \
        -a2 IGH_contig_primary/hifiasm_result/${sample_id}.IGH.a_ctg.fa
    
    #python3 IGLoo_asm.py -rd asm_annotation/hprc/ -id ${sample_id} \
    #    -a1 /home/mlin77/data_blangme2/fasta/hprc-yr1/${sample_id}.1.fa \
    #    -a2 /home/mlin77/data_blangme2/fasta/hprc-yr1/${sample_id}.2.fa
done
