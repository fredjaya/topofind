#!/usr/bin/env python3
import argparse
import os
import subprocess
import re
import glob
import multiprocessing 
from collections import OrderedDict
import sys

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--base_aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-nt', '--num_threads', type=int, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    parser.add_argument('-v', '--verbose', action='store_false', help='Print detailed output')
    parser.add_argument('-r', '--run_nf', action='store_false', help='[DEBUG] Run nextflow processes?')
    parser.add_argument('-s', '--skip_iter', type=int, help='[DEBUG] Only run nextflow processes from n_trees', default=0)
    return parser.parse_args()

def run_nf(run_file, submodel):
    if not args.run_nf or args.skip_iter >= n_trees:
        return
    nf = subprocess.run(["nextflow", "run", f"{repo_path}/main.nf",
        "--base_aln", args.base_aln,
        "--out_dir", args.output_dir,
        "--run_file", run_file, 
        "--submodel", submodel, 
        "--nthreads", str(args.num_threads),
        "-profile", args.executor],
        capture_output=args.verbose)

def get_submodels(n_trees):
    return(",".join(["GTR+FO+G"]*n_trees))

def check_valid_runs(MastResults):
    ''' If there are no new partitions, stop program '''
    none_counter=False
    for key, value in MastResults.items():
        if key.startswith(str(n_trees)):
            if value is not None:
               none_counter=True
    if not none_counter:
        sys.exit("\nStopping program :)\n")

if __name__ == '__main__':
    """
    Initialise variables
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    n_trees = 1 
    bic_improving=True
    
    print(f"\n\
        base_aln:       {args.base_aln}\n\
        output_dir:     {args.output_dir}\n\
        num_threads:    {args.num_threads}\n\
        nf_executor:    {args.executor}\n"
        )

    while bic_improving:
        #check_valid_runs(MastResults)
        n_trees+=1 
        print(f"\nRUNNING PROGRAM WITH {n_trees} TREES\n")
      
        if not os.path.exists(f"{args.output_dir}/NewRuns.tsv"):
            print(f"{args.output_dir}/NewRuns.tsv not found")
            # TODO: return arguments for first iteration and replace print statement
            
        else:
            """
            Read in the updated run names and corresponding partitioned aln 
            """
            with open(f"{args.output_dir}/NewRuns.tsv", "r") as nr:
                prev_runs = nr.read()
            run_nf(f"{args.output_dir}/NewRuns.tsv", get_submodels(n_trees))
            bic_improving = False
