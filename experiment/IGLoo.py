# Wrap up python file for the IGLoo 2nd module
import subprocess
import sys
import os
import argparse

# project modules
from scripts.utils import check_program_install, catch_assert  
from scripts import analyze_recombine
from scripts import summarize_good_report



def align_and_index(ref, input_fasta, prefix):
    command = ' '.join(['minimap2', '-ax', 'map-hifi', ref, input_fasta, '|', 'samtools sort', '>', prefix+'.bam'])
    print(command)
    subprocess.call(command, shell=True)
    command = ' '.join(['samtools', 'index', prefix+'.bam'])
    print(command)
    subprocess.call(command, shell=True)



def main():
    parser = argparse.ArgumentParser(description="The 2nd module (bam/fasta) file analyzer of IGLoo.")
    parser.add_argument('-o', '--out', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-s', '--sample_id', help="Sample ID ['sample'].", default="sample")

    parser.add_argument('-lr', '--list_ref', help='list of reference genome for analysis', nargs='+')
    parser.add_argument('-lb', '--list_bed', help='list of annotated bed files for references', nargs='+')
    parser.add_argument('-f', '--input_fasta', help='input unaligned sequence (.fa/.fq) file for analysis')
    parser.add_argument('-b', '--input_bam', help='input alignment file (.bam) for analysis')
    
    # Module Options
    parser.add_argument('--analyze_read', help='[2] Option to run IGLoo analyze recombination', action='store_true')


    #parser.add_argument('-g', '--genome', help="Path to the reference genome.")
    #parser.add_argument('-v', '--vcf', help="Path to the personal vcf file.")
    #parser.add_argument('-r', '--run_id', help="Run ID ['run'].", default="run")
    ## Process options
    #parser.add_argument('--simulate', help='[1] Option to run biastools simulation.', action='store_true')
    #parser.add_argument('--align',    help='[2] Option to run biastools align.', action='store_true')
    #parser.add_argument('--analyze',  help='[3] Option to run biastools analyze.', action='store_true')
    #parser.add_argument('--predict',  help='[4] Option to predict bias from analysis report.', action='store_true')

    #parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    #parser.add_argument('--force', help="running the program without checking prerequisite programs.", action='store_true')
    ## [1]
    #parser.add_argument('-x', '--coverage', help="Read coverage to simulate [30].", type=int, default=30)
    ## [2]
    #parser.add_argument('-a', '--aligner', help="Aligner to use (bowtie2|bwamem) [bowtie2]", default="bowtie2")
    #parser.add_argument('-b', '--align_index', help="Path to the aligner index (target reference)")
    ## [3]
    #parser.add_argument('-n', '--naive', help= "Option to run the naive assignment method [False].", action='store_true')
    #parser.add_argument('-R', '--real',  help= "Option for performing analysis on real data [False].", action='store_true')
    #parser.add_argument('-d', '--boundary', help= "Boundary to plot the indel balance plot [20]", type=int, default=20)
    #parser.add_argument('-lr', '--list_report', help= "List of bias report to plot the indel balance plot", nargs='+')
    #parser.add_argument('-ld', '--list_run_id', help= "List of run ID for namings in the indel balance plot", nargs='+')
    ## [4]
    #parser.add_argument('-ps', '--sim_report',  help= "Path to the simulation report.")
    #parser.add_argument('-pr', '--real_report', help= "Path to the real read report  [out_dir/sample.real.run.bias].")
    args = parser.parse_args()
    
    ###### Parameters for IGLoo
    path_output = args.out
    sample_id  = args.sample_id
    
    list_ref = args.list_ref
    list_bed = args.list_bed
    input_fasta = args.input_fasta
    input_bam   = args.input_bam

    flag_analyze_read = args.analyze_read
    try:
        assert input_fasta != None or input_bam != None
    except AssertionError:
        catch_assert(parser, "The input BAM or FASTA file should be specified.")

    

    
    #path_ref   = args.genome
    #path_vcf   = args.vcf
    #run_id     = args.run_id
    #
    #flag_simulate = args.simulate
    #flag_align    = args.align
    #flag_analyze  = args.analyze
    #flag_predict  = args.predict

    #path_module = os.path.dirname(__file__) + '/'
    #try:
    #    assert flag_simulate + flag_align + flag_analyze + flag_predict >= 1 
    #except AssertionError:
    #    catch_assert(parser, "At least one of the --simulate/align/analyze/predict option should be specified.")

    #flag_force = args.force
    #thread = args.thread
    #if thread == None:
    #    if sys.platform == "darwin":
    #        result = subprocess.run(["sysctl -n hw.ncpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    #    else:
    #        result = subprocess.run(["nproc"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    #    thread = int(result.stdout.strip())
    #
    #coverage = args.coverage
    #aligner  = args.aligner
    #align_index = args.align_index
    #try:
    #    assert aligner=="bowtie2" or aligner=="bwamem" 
    #except AssertionError:
    #    catch_assert(parser, "Only bowtie2 and bwamem are supported.")

    #flag_naive  = args.naive
    #flag_real   = args.real
    #boundary    = args.boundary
    #list_report = args.list_report
    #list_run_id = args.list_run_id
    #if list_report:
    #    try:
    #        assert len(list_report) == len(list_run_id)
    #    except AssertionError:
    #        catch_assert(parser, "Number of list --list_report and --list_run_id entries are inconsistent.")

    #sim_report  = args.sim_report
    #real_report = args.real_report
    #if flag_predict: 
    #    try:
    #        assert real_report != None
    #    except AssertionError:
    #        catch_assert(parser, "<real_report> should be specified when using --predict")


    #
    ## Checking prerequisite programs are installed
    #if flag_force != True:
    #    check_program_install(["bedtools", \
    #                           "samtools", \
    #                           "bcftools", \
    #                           "bwa", \
    #                           "bowtie2", \
    #                           "gzip", \
    #                           "tabix", \
    #                           "mason_simulator"])

    ## Start running
    #command = "mkdir -p " + path_output
    #subprocess.call(command, shell=True)


    
    print("[IGLoo] Processing " + sample_id + "...")
    subprocess.call("mkdir -p " + path_output, shell=True)
    prefix_out = path_output + '/' + sample_id
    if input_bam != None: # input_bam is specified
        # make sure the input bam file is indexed
        if os.path.isfile(input_bam+'.bai') == False:
            command = ' '.join(['samtools', 'index', input_bam])
            subprocess.call(command, shell=True)
        command = ['-bed', list_bed[0], '-bam', input_bam, '-out', prefix_out + '.1', '-name']
        analyze_recombine.main(command)
        if input_fasta == None:
            #make the IGH fasta_file
            input_fasta = prefix_out + '.IGH.fastq'
            command = ' '.join(['samtools fastq', input_bam, '>', input_fasta])
            subprocess.call(command, shell=True)
    else:
        # Align the raw reads to reference
        align_and_index(list_ref[0], input_fasta, prefix_out+'.0')
        # Analyze the bam file
        command = ['-bed', list_bed[0], '-bam', prefix_out+'.0.bam', '-out', prefix_out + '.1', '-name']
        analyze_recombine.main(command)

    # The loop process
    idx_max=1
    for idx, bed_file in enumerate(list_bed[1:]):
        ref_file = list_ref[idx+1]
        if os.path.isfile(prefix_out+'.'+str(idx+1)+'.bad.rpt'):
            print("[IGLoo] There are bad alignments.")
            command = ' '.join(['sed', "'s/|/ /'", prefix_out+'.'+str(idx+1)+'.bad.rpt', '|', "awk '{print $1}'", '>', \
                                prefix_out+'.'+str(idx+1)+'.bad.reads'])
            print(command)
            subprocess.call(command, shell=True)
            command = ' '.join(['grep -w -A 3 -f',  prefix_out+'.'+str(idx+1)+'.bad.reads', input_fasta, ' --no-group-separator', '>', \
                                prefix_out+'.'+str(idx+1)+'.bad.fastq'])
            print(command)
            subprocess.call(command, shell=True)

            # Align the raw reads to reference
            align_and_index(ref_file, prefix_out+'.'+str(idx+1)+'.bad.fastq', prefix_out+'.'+str(idx+1))
            # Analyze the bam file
            command = ['-bed', bed_file, '-bam', prefix_out+'.'+str(idx+1)+'.bam', '-out', prefix_out + '.' +str(idx+2)]
            analyze_recombine.main(command)
            idx_max = idx+2
        else:
            print("[IGLoo] All remaining alignments are good!")
            break
    list_good_rpt = []
    for idx in range(1, idx_max+1):
        list_good_rpt.append(prefix_out+'.'+str(idx)+'.good.rpt')
    
    fo_summary = open(prefix_out + '.final.rpt', 'w')
    result_good = summarize_good_report.main(['-lr'] + list_good_rpt)
    fo_summary.write('\t'.join([str(ele)for ele in result_good])+'\t')
    
    command=' '.join(['wc -l', prefix_out+'.'+str(idx_max)+'.bad.rpt',  "| awk '{print $1}'"])
    result_bad = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    result_bad = result_bad.stdout.strip()
    fo_summary.write(result_bad+'\n')
    fo_summary.close()




    
if __name__ == "__main__":
    main()


