# Wrap up python file for the IGLoo 2nd module
import subprocess
import sys
import os
import argparse

# project modules
from scripts.utils import check_program_install, catch_assert  
import IGLoo_asm
from scripts import ig_SV_typing
from scripts import merge_personal_ref_2hap
from scripts import merge_haplotypes
from scripts import separate_reads
from scripts import find_upperCase
from scripts import find_coverage_region
from scripts import mask_reassembly



def main():
    parser = argparse.ArgumentParser(description="The 2nd module (bam/fasta) file analyzer of IGLoo.")
    parser.add_argument('-rd', '--result_dir', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-id', '--sample_id',  help="Sample ID ['sample_id'].", default="sample_id")

    parser.add_argument('-fa', '--preprocessed_fasta', help='input preprocessed file', required=True)
    parser.add_argument('-p1', '--parent_1', help='input parent 1 file, should be BAM/FASTA')
    parser.add_argument('-p2', '--parent_2', help='input parent 2 file, should be BAM/FASTA')
    parser.add_argument('-t', '--threads', help='number of threads to use', default=8)
    
    #parser.add_argument('--force', help="running the program without checking prerequisite programs.", action='store_true')
    args = parser.parse_args()
    
    ###### Parameters for IGLoo
    out_dir = args.result_dir
    sample_id  = args.sample_id
    
    input_fasta = args.preprocessed_fasta
    parent_1 = args.parent_1
    parent_2 = args.parent_2
    threads = args.threads

    
    ## Checking prerequisite programs are installed
    #if flag_force != True:
    check_program_install(["chromosome_scaffolder.sh", \
                           "splitScaffoldsAtNs.sh", \
                           "final_polish.sh", \
                           "samtools", \
                           "minimap2", \
                           "jasper.sh", \
                           "jellyfish"])

    # Prepare output directory
    subprocess.call("mkdir -p " + out_dir+'/ref_guide/', shell=True)


    # Run IGLoo_SV_typing
    command = ['-csv1', out_dir+'/asm_annotate/'+sample_id+'.contig_gene.1.csv', \
               '-csv2', out_dir+'/asm_annotate/'+sample_id+'.contig_gene.2.csv', \
               '-f1', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa', \
               '-f2', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa']
    ig_SV_typing.main(command)
    
    # Run merge_personal_ref_2hap
    command = ['-base', 'IGLoo/materials/personalized_ref/referece_base.fa', \
               '-alt',  'IGLoo/materials/personalized_ref/referece_del.fa', \
               '-csv1', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa.rec.csv', \
               '-csv2', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa.rec.csv', \
               '-out1', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.ref.fa', \
               '-out2', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.ref.fa']
    merge_personal_ref_2hap.main(command)


    # Run MaSuRCA
    command = ' '.join(['cp', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.ref.fa', out_dir+'/ref_guide/'])
    subprocess.call(command, shell=True)
    command = ' '.join(['cp', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.ref.fa', out_dir+'/ref_guide/'])
    subprocess.call(command, shell=True)
    
    old_dir = os.getcwd()
    os.chdir(out_dir+'/ref_guide/')

    # Arrange the contigs
    ref_pat = sample_id+'.IGH.asm.hap1.ref.fa'
    ref_mat = sample_id+'.IGH.asm.hap2.ref.fa'
    contigs_1 = sample_id+'.IGH.asm.hap1.fa'
    contigs_2 = sample_id+'.IGH.asm.hap2.fa'

    command = ' '.join(['chromosome_scaffolder.sh', '-r', ref_pat, '-c', '50000', '-i', '99', \
               '-m', '250000', '-q', '../reassembly/'+sample_id+'/'+contigs_1, '-t', str(threads), '-nb'])
    subprocess.call(command, shell=True)
    command = ' '.join(['chromosome_scaffolder.sh', '-r', ref_mat, '-c', '50000', '-i', '99', \
               '-m', '250000', '-q', '../reassembly/'+sample_id+'/'+contigs_2, '-t', str(threads), '-nb'])
    subprocess.call(command, shell=True)

    # split alignments
    arrange_contigs_1 = ref_pat + '.' + contigs_1 + '.split.reconciled.fa'
    arrange_contigs_2 = ref_mat + '.' + contigs_2 + '.split.reconciled.fa'
    command = ' '.join(['splitScaffoldsAtNs.sh', arrange_contigs_1, '1', '>', arrange_contigs_1+'.split'])
    subprocess.call(command, shell=True)
    command = ' '.join(['splitScaffoldsAtNs.sh', arrange_contigs_2, '1', '>', arrange_contigs_2+'.split'])
    subprocess.call(command, shell=True)

    # final polishing
    subprocess.call("mkdir -p " + 'draft_polish_H1', shell=True)
    subprocess.call("mkdir -p " + 'draft_polish_H2', shell=True)
    os.chdir('draft_polish_H1')
    command = ' '.join(['final_polish.sh', '14', '../'+ref_pat, '../'+ref_pat, '../'+arrange_contigs_1])
    subprocess.call(command, shell=True)

    os.chdir('../draft_polish_H2')
    command = ' '.join(['final_polish.sh', '14', '../'+ref_mat, '../'+ref_mat, '../'+arrange_contigs_2])
    subprocess.call(command, shell=True)
    
    # separate reads
    os.chdir('../')

    command = ['-f1', 'draft_polish_H1/14.dir/14.all.polished.fa', '-f2', 'draft_polish_H2/14.dir/14.all.polished.fa', '-out', 'draft_polish.fa']
    merge_haplotypes.main(command)
    command = ' '.join(['minimap2 -ax map-pb draft_polish.fa ' + old_dir+'/'+input_fasta + ' | \
                        samtools sort -o '+sample_id+'.draft.realign.bam'])
    subprocess.call(command, shell=True)
    command = ' '.join(['samtools index', sample_id+'.draft.realign.bam'])
    subprocess.call(command, shell=True)
    command = ['-bam', sample_id+'.draft.realign.bam', '-fasta', old_dir+'/'+input_fasta, '-out', sample_id+'.separate.read']
    separate_reads.main(command)

    
    # Run Jasper for final polishing
    os.chdir('draft_polish_H1')
    command = ' '.join(['jasper.sh', '-t', '16', '-b', '800000000', '-a', '14.dir/14.all.polished.fa', \
                        '-r', '../'+sample_id+'.separate.read.H1.fa', '-k', '25', '-p', '3'])
    subprocess.call(command, shell=True)

    os.chdir('../draft_polish_H2')
    command = ' '.join(['jasper.sh', '-t', '16', '-b', '800000000', '-a', '14.dir/14.all.polished.fa', \
                        '-r', '../'+sample_id+'.separate.read.H2.fa', '-k', '25', '-p', '3'])
    subprocess.call(command, shell=True)
    
    
    # Final masking
    os.chdir('../')

    command = ['-f1', 'draft_polish_H1/14.all.polished.fa.polished.fasta', \
               '-f2', 'draft_polish_H2/14.all.polished.fa.polished.fasta', \
               '-out', 'final_polish.fa']
    merge_haplotypes.main(command)
    command = ' '.join(['minimap2 -ax map-pb final_polish.fa ' + old_dir+'/'+input_fasta + ' | \
              samtools sort -o '+sample_id+'.final.realign.bam'])
    subprocess.call(command, shell=True)
    command = ' '.join(['samtools index', sample_id+'.final.realign.bam'])
    subprocess.call(command, shell=True)

    command = ' '.join(['mkdir -p ../IGLoo_assembly'])
    subprocess.call(command, shell=True)

    command = ['-fasta', 'final_polish.fa', '-out', 'final_polish.upperCase.bed']
    find_upperCase.main(command)
    command = ' '.join(['samtools depth', sample_id+'.final.realign.bam', '-a -g 0x100 -J > '+sample_id+'.realign.rd.log'])
    subprocess.call(command, shell=True)
    command = ['-rd', sample_id+'.realign.rd.log', '-out', sample_id+'.realign.rd.bed']
    find_coverage_region.main(command)
    command = ['-rd', sample_id+'.realign.rd.bed', '-up', 'final_polish.upperCase.bed', '-fa', 'final_polish.fa', '-out', 'final_polish.mask']
    mask_reassembly.main(command)
    command = ' '.join(['cp final_polish.mask.1.fa ../IGLoo_assembly/'+sample_id+'.mask.1.fa'])
    subprocess.call(command, shell=True)
    command = ' '.join(['cp final_polish.mask.2.fa ../IGLoo_assembly/'+sample_id+'.mask.2.fa'])
    subprocess.call(command, shell=True)


    ## Run IGLoo_asm annotation for the draft assembly
    #command = ['-rd', out_dir+'/asm_annotate/', '-id', sample_id, \
    #           '-a1', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa', \
    #           '-a2', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa']
    #IGLoo_asm.main(command)
        
    
    



    
if __name__ == "__main__":
    main()


