# Wrap up python file for the IGLoo 2nd module
import subprocess
import sys
import os
import argparse

# Enable local imports if run directly
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# project modules
try :
    from IGLoo.scripts.utils import check_program_install, catch_assert  
    import IGLoo.IGLoo_asm as IGLoo_asm
except ImportError:
    from scripts.utils import check_program_install, catch_assert  
    import IGLoo_asm




def main(arguments=None):
    parser = argparse.ArgumentParser(description="The 2nd module (bam/fasta) file analyzer of IGLoo.")
    parser.add_argument('-rd', '--result_dir', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-id', '--sample_id',  help="Sample ID ['sample_id'].", default="sample_id")

    parser.add_argument('-fa', '--preprocessed_fasta', help='input preprocessed file', required=True)
    parser.add_argument('-p1', '--parent_bam_1', help='input parent 1 file, should be BAM/FASTA')
    parser.add_argument('-p2', '--parent_bam_2', help='input parent 2 file, should be BAM/FASTA')
    parser.add_argument('-t', '--threads', help='number of threads to use', default=8)
    
    #parser.add_argument('--force', help="running the program without checking prerequisite programs.", action='store_true')
    args = parser.parse_args() if arguments is None else arguments
    
    ###### Parameters for IGLoo
    out_dir = args.result_dir
    sample_id  = args.sample_id
    
    input_fasta = args.preprocessed_fasta
    parent_bam_1 = args.parent_bam_1
    parent_bam_2 = args.parent_bam_2
    threads = args.threads

    
    ## Checking prerequisite programs are installed
    #if flag_force != True:
    check_program_install(["samtools", \
                           "yak", \
                           "hifiasm"])

    # Prepare output directory
    subprocess.call("mkdir -p " + out_dir+'/reassembly/'+sample_id, shell=True)
    subprocess.call("mkdir -p " + out_dir+'/parent/', shell=True)


    
    if parent_bam_1 == None:
        #de novo assembly without parent information
        command = ' '.join(['hifiasm', '-o', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm', '-t '+str(threads), input_fasta])
        subprocess.call(command, shell=True)

        # export to fasta files
        command = ' '.join(['awk', "'/^S/{print \">\"$2;print $3}'", \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.bp.hap1.p_ctg.gfa', '>', \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa'])
        subprocess.call(command, shell=True)
        command = ' '.join(['awk', "'/^S/{print \">\"$2;print $3}'", \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.bp.hap2.p_ctg.gfa', '>', \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa'])
        subprocess.call(command, shell=True)
    else: # Parent information is provided
        # make sure the input files are in the correct format
        if parent_bam_1[-4:] == ".bam":
            command = ' '.join(['samtools', 'fasta', parent_bam_1, '>', parent_bam_1[:-4]+".fa"])
            subprocess.call(command, shell=True)
            parent_bam_1 = parent_bam_1[:-4] + ".fa"
        if parent_bam_2[-4:] == ".bam":
            command = ' '.join(['samtools', 'fasta', parent_bam_2, '>', parent_bam_2[:-4]+".fa"])
            subprocess.call(command, shell=True)
            parent_bam_2 = parent_bam_2[:-4] + ".fa"

        # Trio-binning mode of hifiasm
        command = ' '.join(['yak count -k31 -b37 -t16 -o', out_dir+'/parent/pat.yak', parent_bam_1])
        subprocess.call(command, shell=True)
        command = ' '.join(['yak count -k31 -b37 -t16 -o', out_dir+'/parent/mat.yak', parent_bam_2])
        subprocess.call(command, shell=True)

        command = ' '.join(['hifiasm', '-o', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm', '-t '+str(threads), \
                            '-1', out_dir+'/parent/pat.yak', '-2', out_dir+'/parent/mat.yak', input_fasta])
        subprocess.call(command, shell=True)

        # export to fasta files
        command = ' '.join(['awk', "'/^S/{print \">\"$2;print $3}'", \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.dip.hap1.p_ctg.gfa', '>', \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa'])
        subprocess.call(command, shell=True)
        command = ' '.join(['awk', "'/^S/{print \">\"$2;print $3}'", \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.dip.hap2.p_ctg.gfa', '>', \
                            out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa'])
        subprocess.call(command, shell=True)


    # Run IGLoo_asm annotation for the draft assembly
    #command = ['-rd', out_dir+'/asm_annotate/', '-id', sample_id, \
    #           '-a1', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa', \
    #           '-a2', out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa']
    input_arguments = argparse.Namespace(result_dir=out_dir+'/asm_annotate/', \
                                         sample_id=sample_id, \
                                         assembly_1=out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap1.fa', \
                                         assembly_2=out_dir+'/reassembly/'+sample_id+'/'+sample_id+'.IGH.asm.hap2.fa')
    IGLoo_asm.main(input_arguments)
        
    
    



    
if __name__ == "__main__":
    main()


