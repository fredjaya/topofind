import argparse
import os
import re
import sys

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-f', '--aln_format', choices=['fasta', 'phylip'], help='Alignment format', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-T', '--num_threads', type=str, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    args = parser.parse_args()

if __name__ == '__main__':
    cli()
    ''' 1. Initialise values and dictionary '''
    # A counter for the splitting
    block_name='1'
    
    ''' 2. Run split_aln '''
    run_name="01_split_aln"
    # Run nexflow
    cmd=f"nextflow run main.nf --mode split_aln --aln {args.aln} --aln_format {args.aln_format} --out {args.output_dir} --nthreads {args.num_threads} -profile {args.executor}"
    result=subprocess.run(cmd, capture_output=True)
    # Collect BIC
    # Collect trees
    # Collect partitions
    ''' 3. MAST '''
