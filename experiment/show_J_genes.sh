work_dir="/scratch4/blangme2/maojan/gAIRRsuite/hprc_annotation"
for sample_id in HG002 \
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
do
    echo -n ${sample_id}.1"	"
    python3 show_J_genes.py -bed ${work_dir}/${sample_id}/group_genes.1.bed
    echo -n ${sample_id}.2"	"
    python3 show_J_genes.py -bed ${work_dir}/${sample_id}/group_genes.2.bed
done

