#!/usr/bin/env python3

import argparse
import sys
import os
import logging
import subprocess
from datetime import datetime
import importlib.util

def setup_logging(result_dir):
    """Set up logging configuration"""
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    log_file = os.path.join(result_dir, f"IGLoo_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger('IGLoo')

def import_module(module_name):
    """Dynamically import IGLoo modules"""
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        module_path = os.path.join(script_dir, f"{module_name}.py")
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        raise ImportError(f"Failed to import {module_name}: {str(e)}", exc_info=True)

def run_pipeline_step(step_name, module_name, args, logger):
    """Run a single pipeline step with error handling"""
    logger.info(f"Starting {step_name}")
    try:
        module = import_module(module_name)
        module.main(args)
        logger.info(f"Completed {step_name} successfully")
        return True
    except Exception as e:
        logger.error(f"Error in {step_name}: {str(e)}", exc_info=True)
        return False
    

def main():
    parser = argparse.ArgumentParser(description='IGLoo: Immunoglobulin Locus Analysis Pipeline v0.1.0')
    
    # Main operation mode arguments
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--all', action='store_true', help='Run the complete assembly pipeline in sequence')
    mode_group.add_argument('--asm', action='store_true', help='Run assembly analysis module')
    mode_group.add_argument('--filter', action='store_true', help='Filtering IGH related reads')
    mode_group.add_argument('--read', action='store_true', help='Run read analysis module')
    mode_group.add_argument('--denovo', action='store_true', help='Run de novo assembly module')
    mode_group.add_argument('--refguide', action='store_true', help='Run reference-guided assembly module')

    # Common arguments
    parser.add_argument('-rd', '--result_dir', help="Path to output directory ['result_dir'].", default="result_dir")
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample_id'].", default="sample_id")
    parser.add_argument('-t', '--threads', type=int, default=16, help='Number of threads')

    # Module-specific arguments
    # Assembly module
    parser.add_argument('-a1', '--assembly_1', help='input assembly haplotype 1 file (.fa) for analysis')#, required=True)
    parser.add_argument('-a2', '--assembly_2', help='input assembly haplotype 2 file (.fa) for analysis')

    # Filtering module
    parser.add_argument('-ref', '--ref_genome', help='reference genome for filtering [GRCh38/CHM13]', default="GRCh38")
    parser.add_argument('--unmap', action='store_true', help='include unmapped reads in the filterred bam file')

    # Read module
    parser.add_argument('-lr', '--list_ref', help='list of reference genome for analysis', nargs='+')
    parser.add_argument('-lb', '--list_bed', help='list of annotated bed files for references', nargs='+')
    parser.add_argument('-f', '--input_fasta', help='input unaligned sequence (.fa/.fq) file for analysis')
    parser.add_argument('-b', '--input_bam', help='input alignment file (.bam) for analysis')
    parser.add_argument('-ont', '--nanopore', action='store_true', help='flag for nanopore data')

    # De novo & Reference-guided assembly module
    parser.add_argument('-fa', '--preprocessed_fasta', help='input preprocessed file')#, required=True) # should be --read module output
    parser.add_argument('-pb1', '--parent_bam_1', help='input parent 1 file, should be BAM/FASTA')
    parser.add_argument('-pb2', '--parent_bam_2', help='input parent 2 file, should be BAM/FASTA')
    parser.add_argument('-py1', '--parent_yak_1', help='input parent 1 file, should be yak format [result_dir/parent/pat.yak].')
    parser.add_argument('-py2', '--parent_yak_2', help='input parent 2 file, should be yak format [result_dir/parent/mat.yak].')
    
    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.result_dir)
    logger.info("Starting IGLoo pipeline")

    #for element in args:
    #    print(element)

    try:
        success = False
        error_msg = None
        if args.all:
            # Run full pipeline in sequence
            steps = [
                ("Read Analysis", "IGLoo_read"),
                ("De novo Assembly", "IGLoo_ReAsm"),
                ("Reference-guided Assembly", "IGLoo_ReAsm2"),
            ]
            
            for step_name, module_name in steps:
                if not run_pipeline_step(step_name, module_name, args, logger):
                    error_msg = f"Failed in {step_name} step"
                    break
            else:  # This runs if no break occurred
                success = True
        else:
            # Run individual modules
            if args.filter:
                assert args.input_bam is not None, "input bam file is required for filtering"
                assert args.ref_genome in ["GRCh38", "CHM13"], "reference genome should be either GRCh38 or CHM13"
                try:
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    if args.ref_genome == "GRCh38":
                        module_path = os.path.join(script_dir, f"scripts/collect_IGH_from_grch38.sh")
                    else:
                        module_path = os.path.join(script_dir, f"scripts/collect_IGH_from_chm13.sh")
                    str_unmap = "--unmap" if args.unmap else ""
                    subprocess.run(' '.join(['bash', module_path, args.sample_id, args.input_bam, args.result_dir, str_unmap]), shell=True)
                    success = True
                except Exception as e:
                    logger.error(f"Error in Filtering: {str(e)}", exc_info=True)
                    error_msg = f"Filtering failed: {str(e)}"
            elif args.read:
                success = run_pipeline_step("Read Analysis", "IGLoo_read", args, logger)
                if not success: error_msg = "Read Analysis failed"
            elif args.asm:
                success = run_pipeline_step("Assembly Analysis", "IGLoo_asm", args, logger)
                if not success: error_msg = "Assembly Analysis failed"
            elif args.denovo:
                success = run_pipeline_step("De novo Assembly", "IGLoo_ReAsm", args, logger)
                if not success: error_msg = "De novo Assembly failed"
            elif args.refguide:
                success = run_pipeline_step("Reference-guided Assembly", "IGLoo_ReAsm2", args, logger)
                if not success: error_msg = "Reference-guided Assembly failed"

        if success:
            logger.info("IGLoo pipeline completed successfully")
        else:
            raise Exception(error_msg if error_msg else "Pipeline failed with unexpected error")

    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
