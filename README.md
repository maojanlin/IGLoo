
_Updated: Jun 08, 2024_
# IGLoo
Analyzing the Immunoglobulin (IG) HiFi read data and assemblies derived from Lymphoblastoid cell lines (LCLs).



## Usage

### Profiling Assemblies
```
$ python3 IGLoo/IGLoo_asm.py [-h] [-rd RESULT_DIR] [-id SAMPLE_ID] -a1 ASSEMBLY_1 [-a2 ASSEMBLY_2]
```


## Examples
Running the ```IGLoo --asm```
```
$ python3 IGLoo/IGLoo_asm.py -rd example/asm_out/ -id HG005 -a1 example/HG005.hprc.asm.1.IGH.fa -a2 example/HG005.hprc.asm.2.IGH.fa
```
