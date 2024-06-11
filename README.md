
_Updated: Jun 10, 2024_
# IGLoo
Analyzing the Immunoglobulin (IG) HiFi read data and assemblies derived from Lymphoblastoid cell lines (LCLs).



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
