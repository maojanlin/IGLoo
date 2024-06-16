#python3 IGLoo/IGLoo_ReAsm.py -rd example/ReAsm_out/ -id HG005 \
#                             -fa example/read_out2/processed_fasta/HG005.split.enrich.fa \
#                             -p1 example/HG006.final.IGH.bam \
#                             -p2 example/HG007.final.IGH.bam

export PYTHONPATH=/home/mlin77/lib/python3.9
python3 IGLoo/IGLoo_ReAsm2.py -rd example/ReAsm_out/ -id HG005 \
                              -fa example/read_out2/processed_fasta/HG005.split.enrich.fa \
                              -p1 example/HG006.final.IGH.bam \
                              -p2 example/HG007.final.IGH.bam
