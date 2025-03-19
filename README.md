
_Updated: Mar 28, 2025_

[![DOI](https://zenodo.org/badge/631047974.svg)](https://doi.org/10.5281/zenodo.15048412)

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
### Python packages
- numpy
- pysam
- pandas
- matplotlib
- seaborn

## Installation
- [pip](https://pypi.org/project/bio-IGLoo/)
```
pip install bio-IGLoo
```
- [Github](https://github.com/maojanlin/IGLoo.git)
```
git clone https://github.com/maojanlin/IGLoo.git
cd IGLoo
```
Though optional, it is a good practice to install a virtual environment to manage the dependancies:

```
python -m venv venv
source venv/bin/activate
```
Now a virtual environment (named venv) is activated. Install biastools:

```
python setup.py install
```


## Usage

### Fetch IGH related reads from a bam file
User can use the utility scripts ```IGLoo/scripts/collect_IGH_from_grch38.sh``` or ```IGLoo/scripts/collect_IGH_from_chm13.sh``` to subset the IGH related alignments from full WGS alignments.
```
$ bash IGLoo/scripts/collect_IGH_from_grch38.sh HG005 ./HG005_aligned_GRCh38_winnowmap.sorted.bam ./example/
```
Or directly from IGLoo
```
$ IGLoo --filter [-rd RESULT_DIR] [-id SAMPLE_ID] -b INPUT_BAM --ref_genome <"GRCh38"/"CHM13">
```


### Profiling Assemblies
```
$ IGLoo --asm [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -a1 ASSEMBLY_1 [-a2 ASSEMBLY_2]
```

### Profiling HiFi read data
```
$ IGLoo --read [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] \
               -b  INPUT_BAM \
               -f  INPUT_FASTA \
               -lb BED_1 [BED_2 BED_3 ... ] \
               -lr REF_1 [REF_2 REF_3 ... ]
```

Note that at least one of the ```INPUT_BAM``` and ```INPUT_FASTA``` should be specified.  The reference used in the ```INPUT_BAM``` should be put in the first position.

### Reassemble personal assemblies in the IGH

```
$ IGLoo --denovo   [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -fa PREPROCESSED_FASTA [-pb1 PARENT_1_BAM] [-pb2 PARENT_2_BAM] [-t THREADS]
$ IGLoo --refguide [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -fa PREPROCESSED_FASTA [-py1 PARENT_1_YAK] [-py2 PARENT_2_YAK] [-t THREADS]
```

The ```--denovo``` command assembles the draft assembly with hifiasm.  The ```--refguide``` types the SV from draft assembly to generate personal references.  MaSuRCA was then applied to generate the final assembly.
Note that to successfully runs the JASPER and Jellyfish in MaSuRCA, ```$PYTHONPATH``` needs to be set with ```export PYTHONPATH=/home/$USER/lib/python3.9```.


## Examples
### Running the ```IGLoo --asm```
```
$ python3 IGLoo/IGLoo_asm.py -rd example/asm_out/ -id HG005 -a1 example/HG005.hprc.asm.1.IGH.fa -a2 example/HG005.hprc.asm.2.IGH.fa
```

### Running the ```IGLoo --read```


Then perform the IGLoo --read analysis on the subset bam file ```HG005.hprc.IGH.bam```.

```
$ IGLoo --read -id HG005 \
               -rd example/read_out/ \
               -b  example/HG005.hprc.IGH.bam \
               -lb IGLoo/materials/gene_annotations/GRCh38/grch38_IGH.bed \
                   IGLoo/materials/gene_annotations/hg19_IGH.bed \
                   IGLoo/materials/gene_annotations/chm13_IGH.bed \
               -lr path_to_grch38.fa \
                   path_to_grch37.fa \
                   path_to_chm13.fa
```

The final results will be generated in ```example/read_out/pc_report/```, including ```HG005.split.rpt```, which summarizes overall recombination events and their frequencies (read counts);  ```HG005.split.detail.rpt```, which lists each event alongside its corresponding read name; and ```HG005.pie_chart.pdf```, which shows a pie chart of all recombination events.


#### Running with nanopore data
Use the flag ```--nanopore``` for analyzing the recombination events using nanopore sequence data.
```
$ IGLoo --read --nanopore \
               -id HG005 \
               ...
```



### Running the ```IGLoo --ReAsm```
```
IGLoo --denovo -rd example/ReAsm_out/ -id HG005 \
               -fa example/read_out2/processed_fasta/HG005.split.enrich.fa \
               -pb1 example/HG006.final.IGH.bam \
               -pb2 example/HG007.final.IGH.bam

export PYTHONPATH=/home/$USER/lib/python3.9
IGLoo --refguide -rd example/ReAsm_out/ -id HG005 \
                 -fa example/read_out2/processed_fasta/HG005.split.fa
```

Note that $PYTHONPATH needs to be specified to run JASPER and Jellyfish.

