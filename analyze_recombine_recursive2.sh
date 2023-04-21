echo "sample_id	#CVDJ	#VDJ	#total_read	#max_read	#unrecomb	#bad_reads" > final_summary.2.rpt

summarize() {
    echo "[AIRRsembly] Processing ${1}..."
    python3 analyze_recombine.py -bed materials/gene_annotations/GRCh38/grch38_IGH.bed \
                                 -bam extract_hprc_bam/${1}.hprc.IGH.bam \
                                 -out recombine_results/${1}.1
    
    if [ -s recombine_results/${1}.1.bad.rpt ]; then
        # The file is not-empty.
        echo "[AIRRsembly] There are bad alignments."
        sed 's/|/ /' recombine_results/${1}.1.bad.rpt | awk '{print $1}' > \
                     recombine_results/${1}.1.bad.names
        grep -w -A 3 -f recombine_results/${1}.1.bad.names sequence_files/${1}.hprc.IGH.fastq \
                        --no-group-separator > recombine_results/${1}.1.bad.fastq
        
        hg19="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa"
        hg19_bed="materials/gene_annotations/hg19_IGH.bed"
        echo "[AIRRsembly] Realignment on hg19..."
        #minimap2 -ax map-hifi ${hg19} recombine_results/${1}.1.bad.fastq > \
        #                              recombine_results/${1}.1.bad.bam
        #samtools sort recombine_results/${1}.1.bad.bam > recombine_results/${1}.1.bad.sorted.bam
        #samtools index recombine_results/${1}.1.bad.sorted.bam
        
        # second round analysis
        python3 analyze_recombine.py -bed ${hg19_bed} \
                                     -bam recombine_results/${1}.1.bad.sorted.bam \
                                     -out recombine_results/${1}.2
        if [ -s recombine_results/${1}.2.bad.rpt ]; then
            # The file is not-empty.
            echo "[AIRRsembly] There are bad alignments."
            sed 's/|/ /' recombine_results/${1}.2.bad.rpt | awk '{print $1}' > \
                         recombine_results/${1}.2.bad.names
            grep -w -A 3 -f recombine_results/${1}.2.bad.names sequence_files/${1}.hprc.IGH.fastq \
                            --no-group-separator > recombine_results/${1}.2.bad.fastq
            
            hg01258_ref="/scratch4/blangme2/maojan/AIRRssembly/HG01258_unique/grch37.HG01258.unique.fa"
            hg01258_ref_bed="/scratch4/blangme2/maojan/AIRRssembly/HG01258_unique/grch37.HG01258.unique.IGH.bed"
            minimap2 -ax map-hifi ${hg01258_ref} recombine_results/${1}.2.bad.fastq > \
                                                 recombine_results/${1}.2.bad.bam
            samtools sort recombine_results/${1}.2.bad.bam > recombine_results/${1}.2.bad.sorted.bam
            samtools index recombine_results/${1}.2.bad.sorted.bam
            
            # third round analysis
            python3 analyze_recombine.py -bed ${hg01258_ref_bed} \
                                         -bam recombine_results/${1}.2.bad.sorted.bam \
                                         -out recombine_results/${1}.3
            
            a=$(python3 summarize_good_report.py -lr recombine_results/${1}.1.good.rpt \
                                                     recombine_results/${1}.2.good.rpt \
                                                     recombine_results/${1}.3.good.rpt)
            b=$(wc -l recombine_results/${1}.3.bad.rpt | awk '{print $1}')
            echo ${1}"	""${a}""	"${b} >> final_summary.2.rpt
        else
            a=$(python3 summarize_good_report.py -lr recombine_results/${1}.1.good.rpt \
                                                     recombine_results/${1}.2.good.rpt)
            b=$(wc -l recombine_results/${1}.2.bad.rpt | awk '{print $1}')
            echo ${1}"	""${a}""	"${b} >> final_summary.2.rpt
        fi
    else
        # The file is empty.
        echo "[AIRRsembly] All alignments are good!"
        a=$(python3 summarize_good_report.py -lr recombine_results/${1}.1.good.rpt)
        echo ${1}"	""${a}""	0"    >> final_summary.2.rpt
    fi
}
export -f summarize
parallel summarize :::  HG002 \
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
