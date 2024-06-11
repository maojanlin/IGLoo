# Wrap up python file for the IGLoo --read module
import subprocess
import sys
import os
import argparse

# project modules
from scripts.utils import check_program_install, catch_assert  
from scripts import analyze_pacbio_refs
from scripts import enrich_DJ_read



def align_and_index(ref, input_fasta, prefix):
    command = ' '.join(['minimap2', '-ax', 'map-hifi', ref, input_fasta, '|', 'samtools sort', '>', prefix+'.bam'])
    print(command)
    subprocess.call(command, shell=True)
    command = ' '.join(['samtools', 'index', prefix+'.bam'])
    print(command)
    subprocess.call(command, shell=True)



def main():
    parser = argparse.ArgumentParser(description="The 2nd module (bam/fasta) file analyzer of IGLoo.")
    parser.add_argument('-rd', '--result_dir', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample_id'].", default="sample_id")

    parser.add_argument('-lr', '--list_ref', help='list of reference genome for analysis', nargs='+')
    parser.add_argument('-lb', '--list_bed', help='list of annotated bed files for references', nargs='+')
    parser.add_argument('-f', '--input_fasta', help='input unaligned sequence (.fa/.fq) file for analysis')
    parser.add_argument('-b', '--input_bam', help='input alignment file (.bam) for analysis')
    args = parser.parse_args()
    
    ###### Parameters for IGLoo
    path_output = args.result_dir
    sample_id  = args.sample_id
    
    list_ref = args.list_ref
    list_bed = args.list_bed
    input_fasta = args.input_fasta
    input_bam   = args.input_bam

    try:
        assert input_fasta != None or input_bam != None
    except AssertionError:
        catch_assert(parser, "The input BAM or FASTA file should be specified.")

    ## Checking prerequisite programs are installed
    #if flag_force != True:
    check_program_install(["samtools", \
                           "minimap2"])

    ## Start running
    print("[IGLoo --read] Processing " + sample_id + "...")
    subprocess.call("mkdir -p " + path_output, shell=True)
    prefix_out = path_output + '/' + sample_id
    if input_bam == None: # Input is a fasta file
        align_and_index(list_ref[0], input_fasta, prefix_out+'.0')
        input_bam = prefix_out + '.0.bam'
    if os.path.isfile(input_bam+'.bai') == False: # make sure there is a bam file
        command = ' '.join(['samtools', 'index', input_bam])
        subprocess.call(command, shell=True)

    # make fasta and fastq files
    command = ' '.join(['samtools fastq', input_bam, '>', prefix_out + '.0.fastq'])
    subprocess.call(command, shell=True)
    command = ' '.join(['samtools fasta', input_bam, '>', prefix_out + '.0.fasta'])
    subprocess.call(command, shell=True)

    # Alignment process
    for idx, ref in enumerate(list_ref):
        align_and_index(ref, prefix_out + '.0.fastq', prefix_out+'.'+str(idx+1))

    # make directory for processed fasta files
    command = 'mkdir -p ' + path_output + '/processed_fasta'
    subprocess.call(command, shell=True)
    command = 'mkdir -p ' + path_output + '/pc_report'
    subprocess.call(command, shell=True)

    # run analyze_pacbio_refs
    processed_bam = [input_bam] + [prefix_out+'.'+str(idx+1)+'.bam' for idx in range(len(list_ref))]
    processed_bed = [list_bed[0]] + list_bed

    command = ['-lbed'] + processed_bed + ['-lbam'] + processed_bam + \
              ['-fasta',  prefix_out+'.0.fasta', \
               '-out_fa', path_output+'/processed_fasta/'+sample_id+'.split.fa', \
               '-out_dt', path_output+'/pc_report/'+sample_id+'.split.detail.rpt', \
               '-out',    path_output+'/pc_report/'+sample_id+'.split.rpt']
    analyze_pacbio_refs.main(command)

    # run read preprocessing
    command = ['-fasta', path_output+'/processed_fasta/'+sample_id+'.split.fa', \
               '-out',   path_output+'/processed_fasta/'+sample_id+'.split.enrich.fa']
    enrich_DJ_read.main(command)






    
if __name__ == "__main__":
    main()



