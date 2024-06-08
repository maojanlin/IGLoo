# Wrap up python file for the IGLoo 1st module
import subprocess
import sys
import os
import argparse

# project modules
from scripts.utils import check_program_install, catch_assert  
from scripts import contig_gene_table
from scripts import count_gene_and_allele
from scripts import plot_asm



def main():
    parser = argparse.ArgumentParser(description="The 1st module (assembly) analyzer of IGLoo, only IGH is supported now.")
    parser.add_argument('-rd', '--result_dir', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample_id'].", default="sample_id")

    parser.add_argument('-a1', '--assembly_1', help='input assembly H1 file (.fa) for analysis', required=True)
    parser.add_argument('-a2', '--assembly_2', help='input assembly H2 file (.fa) for analysis')
    args = parser.parse_args()
    

    ###### Parameters for IGLoo
    out_dir = args.result_dir
    sample_id  = args.sample_id
    
    assembly_1 = args.assembly_1
    assembly_2 = args.assembly_2

    # set the environment
    path_module   = os.path.dirname(__file__)
    path_material = os.path.dirname(__file__) + '/materials/'
    
    work_dir = out_dir + '/annotate/'
    subprocess.call("mkdir -p " + out_dir, shell=True)
    subprocess.call("mkdir -p " + work_dir, shell=True)
    
    command = ["gAIRR_annotate", "-wd", work_dir, "-id", sample_id, "-lc", "IG", "-a1", assembly_1] 
    if assembly_2:
        command += ["-a2", assembly_2]
    subprocess.call(command)
    
    # Summarize the IG genes
    command = ['-bed1', work_dir+sample_id+'/group_genes.1.bed', '-out1', out_dir+'/'+sample_id+'.contig_gene.1.csv', \
               '-bed2', work_dir+sample_id+'/group_genes.2.bed', '-out2', out_dir+'/'+sample_id+'.contig_gene.2.csv', \
               '-target', path_material+"/IGH_functional.txt", '--summary_num', out_dir+'/summary.txt']
    contig_gene_table.main(command)

    # Plot the gene count
    command = ['-csv', out_dir+'/'+sample_id+'.contig_gene.1.csv', '-id', sample_id + '.1', '-out', out_dir+'/'+sample_id+'.contig_gene.1.pdf']
    plot_asm.main(command)
    command = ['-csv', out_dir+'/'+sample_id+'.contig_gene.2.csv', '-id', sample_id + '.2', '-out', out_dir+'/'+sample_id+'.contig_gene.2.pdf']
    plot_asm.main(command)
    



    
if __name__ == "__main__":
    main()


