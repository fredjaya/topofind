import argparse
import os
import subprocess
import re
import glob
import multiprocessing 
from Bio import Phylo

class IterationResults:
    """
    Class to save the output of a split_aln and MAST iteration
    """
    def __init__(self):
        self.in_trees = None 
        self.mast_bic = None
        self.mast_alns = None 

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-f', '--aln_format', choices=['fasta', 'phylip'], help='Alignment format', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-T', '--num_threads', type=int, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    parser.add_argument('-v', '--verbose', action='store_false', help='Print detailed output')
    args = parser.parse_args()
    return(args)

def run_nf(aln, run_name):
    subprocess.run(["nextflow", "run", f"{repo_path}/python.nf", 
        "--aln", args.aln, 
        "--aln_format", args.aln_format, 
        "--run_name", run_name, 
        "--out", args.output_dir, 
        "--nthreads", str(args.num_threads),
        "-profile", args.executor], capture_output=args.verbose)

def pool_iterations(alns):
    with multiprocessing.Pool() as pool:
        pool.map(run_iteration, alns)

def just_file(path):
    base=os.path.basename(path)
    return(os.path.splitext(base)[0])

def get_bic(iqtree_out_path):
    line=subprocess.run(["grep", "BIC", iqtree_out_path], capture_output=True, text=True)
    line=line.stdout
    bic=re.sub("^.+: ", "", line)
    return(float(bic))

def parse_outputs(aln, run_name):
    aln_name=just_file(aln)
    out_path=f"{args.output_dir}/{aln_name}/{run_name}/"
   
    # Collect trees used as input to mast
    in_trees=Phylo.parse(f"{out_path}/*.treefile", "newick")

    # Record MAST BIC
    iqtree_out_path=f"{out_path}/t2_mast_tr.iqtree"
    bic=get_bic(iqtree_out_path)

    # Collect fastas post-HMM assignment
    alns=glob.glob(f"{out_path}/*.fas")

    return(in_trees, bic, alns)

def run_iteration(aln, run_name, Results):
    """
    Run split_aln and MAST == one iteration
    """
    run_nf(aln, run_name)
    
    """
    Parse and save the output of a run to a class
    """
    iter_results = IterationResults()
    iter_results.bic, iter_results.alns = parse_outputs(aln, run_name)
    
    """
    Store in the main Results dictionary
    """
    Results[run_name] = iter_results
    return(Results)

if __name__ == '__main__':
    """
    Initialise variables and Results dictionary
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    aln=args.aln
    n_iter = 1
    run_name = f"iter{n_iter}"
    Results = {}

    """
    Take an input alignment, construct two subtrees from +R2 site assignment,
    input subtrees to MAST
    """
    print("\nPartitioning sites according to rate category assignments, and making subtrees...\n")
    Results = run_iteration(aln, run_name, Results)
    print(vars(Results[run_name]))

    n_iter+=1
    run_name=f"iter{n_iter}"
    print(f"\nCommencing iteration {n_iter}")
