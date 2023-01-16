import argparse
import os
import subprocess
import re
import glob
import multiprocessing 

class IterationResults:
    """
    Class to save the output of a split_aln and MAST iteration
    """
    def __init__(self):
        self.bic = None
        self.trees = None 
        self.partitions = None
        self.alns = None 

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-f', '--aln_format', choices=['fasta', 'phylip'], help='Alignment format', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-T', '--num_threads', type=str, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    args = parser.parse_args()
    return(args)

def run_iteration(aln, run_name):
    cmd=f"nextflow run {repo_path}/python.nf --aln {args.aln} --aln_format {args.aln_format} --run_name --out {args.output_dir} --nthreads {args.num_threads} -profile {args.executor}"
    print("\nPartitioning sites according to rate category assignments, and making subtrees...")
    print(f"\n{cmd}\n")
    os.system(cmd)

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

def parse_outputs(aln):
    aln_name=just_file(aln)
    out_path=f"{args.output_dir}/{aln_name}/02_mast/"
    
    # Record BIC
    iqtree_out_path=f"{out_path}/t2_mast_tr.iqtree"
    bic=get_bic(iqtree_out_path)

    # Collect fastas
    alns=glob.glob(f"{args.output_dir}/{aln_name}/02_mast/*.fas")

    # Collect trees
    # Collect partitions

    return(bic, alns)

if __name__ == '__main__':
    ''' 1. Initialise values and dictionary '''
    args=set_args()
    repo_path=os.path.dirname(__file__)
    run_name='1A'
    aln=args.aln
    Results = {}

    if run_name == "1A":
        """
        Run split_aln and MAST == one iteration
        """
        run_iteration(aln, run_name)

        """
        Parse and save the output of a run to a class
        """
        iter_results = IterationResults()
        iter_results.bic, iter_results.alns = parse_outputs(aln)

        """
        Store in the main Results dictionary
        """
        Results[run_name] = iter_results
        print(vars(Results["1A"]))

    else:
        pool_iterations(aln, run_name)
