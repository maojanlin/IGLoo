pacbio_reference_clip_preprocess() {
    # parameters
    experiment_ID="trio_pacbio_3rd_2"
    out_dir="/home/mlin77/vast/maojan/IGLoo_reassembly/${experiment_ID}/"
    mkdir -p ${out_dir}
    mkdir -p ${out_dir}/processed_fasta/
    mkdir -p ${out_dir}/reassembly/
    mkdir -p ${out_dir}/asm_annotate/
    
    # more parameters
    sample_ID=${1}
    sample_fasta="depleted_fasta/${sample_ID}.hprc.IGH.fasta"
    sample_bam_0="extract_hprc_bam/${sample_ID}.hprc.IGH.bam"
    sample_bam_1="IGH_mappings_winnowmap/${sample_ID}.IGH.grch38.winnowmap.bam"
    sample_bam_2="IGH_mappings_winnowmap/${sample_ID}.IGH.grch37.winnowmap.bam"
    sample_bam_3="IGH_mappings_winnowmap/${sample_ID}.IGH.chm13.winnowmap.bam"

    parent_1_yak="parent_data_clipped/${sample_ID}/pat.yak"
    parent_2_yak="parent_data_clipped/${sample_ID}/mat.yak"

    bed_1="materials/gene_annotations/GRCh38/grch38_IGH.bed"
    bed_2="materials/gene_annotations/hg19_IGH.bed"
    bed_3="materials/gene_annotations/chm13_IGH.bed"
    
    echo "[Deplete ReAssembly] Processing ${sample_ID}..."
    python3 analyze_pacbio_refs.py -lbed ${bed_1} \
                                         ${bed_1} \
                                         ${bed_2} \
                                         ${bed_3} \
                                   -lbam ${sample_bam_0} \
                                         ${sample_bam_1} \
                                         ${sample_bam_2} \
                                         ${sample_bam_3} \
                                   -fasta ${sample_fasta} \
                                   -out_fa ${out_dir}/processed_fasta/${sample_ID}.split.fa \
                                   -out_dt pc_report/${sample_ID}.split.detail.rpt \
                                   -out    pc_report/${sample_ID}.split.rpt
    
    #python3 enrich_DJ_read.py -fasta ${out_dir}/processed_fasta/${sample_ID}.split.fa \
    #                          -out   ${out_dir}/processed_fasta/${sample_ID}.split.enrich.fa
    #
    ## hifiasm
    #echo "[AIRRsembly] Assemble with hifiasm trio-binning method..."
    #mkdir -p ${out_dir}/reassembly/${sample_ID}
    #/home/mlin77/scr4_blangme2/maojan/hifiasm/hifiasm -o ${out_dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm -t 8 \
    #                                                  -1 ${parent_1_yak} \
    #                                                  -2 ${parent_2_yak} \
    #                                                  ${out_dir}/processed_fasta/${sample_ID}.split.enrich.fa
    #
    #awk '/^S/{print ">"$2;print $3}' ${out_dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.dip.hap1.p_ctg.gfa > ${out_dir}/reassembly/${sample_ID}.IGH.asm.hap1.fa
    #awk '/^S/{print ">"$2;print $3}' ${out_dir}/reassembly/${sample_ID}/${sample_ID}.IGH.asm.dip.hap2.p_ctg.gfa > ${out_dir}/reassembly/${sample_ID}.IGH.asm.hap2.fa
    #
    ## IGLoo --genotype
    #echo "[AIRRsembly] IGLoo --genotype..."
    #python3 IGLoo_asm.py -rd ${out_dir}/asm_annotate/ -id ${sample_ID} \
    #    -a1 ${out_dir}/reassembly/${sample_ID}.IGH.asm.hap1.fa \
    #    -a2 ${out_dir}/reassembly/${sample_ID}.IGH.asm.hap2.fa
}
export -f pacbio_reference_clip_preprocess 
parallel pacbio_reference_clip_preprocess ::: \
    HG002   \
    HG00438 \
    HG005   \
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
    NA20129 \
    NA21309

