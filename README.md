
_Updated: Jun 16, 2024_
# IGLoo
Analyzing the Immunoglobulin (IG) HiFi read data and assemblies derived from Lymphoblastoid cell lines (LCLs).


## Prerequisite programs
- samtools=v1.11
- minimap2=v2.24
- yak=v0.1-r56
- hifiasm=v0.18.2-r467
- gAIRR-suite=v0.2.0
- MaSuRCA=v4.1.0
- JASPER=v1.0.2


## Usage

### Profiling Assemblies
```
$ python3 IGLoo/IGLoo_asm.py [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -a1 ASSEMBLY_1 [-a2 ASSEMBLY_2]
```

### Profiling HiFi read data
```
$ python3 IGLoo/IGLoo_read.py [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] \
                              -b  INPUT_BAM \
                              -f  INPUT_FASTA \
                              -lb BED_1 [BED_2 BED_3 ... ] \
                              -lr REF_1 [REF_2 REF_3 ... ]
```

Note that at least one of the ```INPUT_BAM``` and ```INPUT_FASTA``` should be specified.  The reference used in the ```INPUT_BAM``` should be put in the first position.

### Reassemble personal assemblies in the IGH

```
$ python3 IGLoo/IGLoo_ReAsm.py  [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -fa PREPROCESSED_FASTA [-p1 PARENT_1] [-p2 PARENT_2] [-t THREADS]
$ python3 IGLoo/IGLoo_ReAsm2.py [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -fa PREPROCESSED_FASTA [-t THREADS]
```

The ```IGLoo_ReAsm.py``` assembles the draft assembly with hifiasm.  The ```IGLoo_ReAsm2.py``` types the SV from draft assembly to generate personal references.  MaSuRCA was then applied to generate the final assembly.



## Examples
Running the ```IGLoo --asm```
```
$ python3 IGLoo/IGLoo_asm.py -rd example/asm_out/ -id HG005 -a1 example/HG005.hprc.asm.1.IGH.fa -a2 example/HG005.hprc.asm.2.IGH.fa
```

Running the ```IGLoo --read```
```
$ python3 IGLoo/IGLoo_read.py -id HG005 \
                              -rd example/read_out/ \
                              -b  example/HG005.hprc.IGH.bam \
                              -lb IGLoo/materials/gene_annotations/GRCh38/grch38_IGH.bed \
                                  IGLoo/materials/gene_annotations/hg19_IGH.bed \
                                  IGLoo/materials/gene_annotations/chm13_IGH.bed \
                              -lr path_to_grch38.fa \
                                  path_to_grch37.fa \
                                  path_to_chm13.fa
```

Running the ```IGLoo --ReAsm```
```
python3 IGLoo/IGLoo_ReAsm.py -rd example/ReAsm_out/ -id HG005 \
                             -fa example/read_out2/processed_fasta/HG005.split.enrich.fa \
                             -p1 example/HG006.final.IGH.bam \
                             -p2 example/HG007.final.IGH.bam

export PYTHONPATH=/home/user/lib/python3.9
python3 IGLoo/IGLoo_ReAsm2.py -rd example/ReAsm_out/ -id HG005 \
                              -fa example/read_out2/processed_fasta/HG005.split.enrich.fa
```

Note that $PYTHONPATH needs to be specified to run JASPER and Jellyfish.

