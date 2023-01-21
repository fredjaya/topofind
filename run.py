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
        "-profile", args.executor],
        capture_output=args.verbose)

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

def collect_subtrees(tree_file, tree_names):
    # Possibly need to read in partitions for the trees as well
    with open (tree_file, "r") as file:
        subtrees = {}
        for tree_name, line in zip(tree_names, file):
            subtrees[tree_name] = line.strip()
    return(subtrees)

def collect_alns(aln_glob, tree_names):
    # Add primes to end of tree_names as partitions are updated after MAST
    tree_primes = [s + "′" for s in tree_names]
    aln_files={}
    for aln_file, tree_name in zip(sorted(aln_glob), tree_primes):
        aln_files[tree_name] = aln_file
    return(aln_files)

def recurse_trees(tree_list):
    """
    This operates on a list of existing trees and makes a new
    list of trees for the next iteration
    """
    if tree_list == None:
        return(['A', 'B'])
    else:
        t0 = str(new_list[0])
        t1 = f"{t0}A" 
        t2 = f"{t0}B"
        new_list += [t1, t2]
        new_list.remove(t0)
        return(new_list)

def parse_outputs(aln, run, new_trees):
    aln_name=just_file(aln)
    out_path=f"{args.output_dir}/{aln_name}/{run}"

    # Record MAST BIC
    iqtree_out_path=f"{out_path}/t2_mast_tr.iqtree"
    bic=get_bic(iqtree_out_path)

    # Collect subtrees used as input to mast
    subtrees=collect_subtrees(f"{out_path}/concat.treefile", new_trees)

    # Collect fastas post-HMM assignment
    ## Get partitions
    ## Get alignment files
    alns=collect_alns(glob.glob(f"{out_path}/*.fas"), new_trees)
    return(bic, subtrees, alns)

def run_iteration(aln, tree_names, n_iter, Results):
    """
    Update names for current iteration
    """
    n_iter+=1
    new_trees = recurse_trees(tree_names)
    run=f"{n_iter}_{'_'.join(new_trees)}"

    print(f"\nRunning {run}...\n")
    """
    Run split_aln and MAST == one iteration
    """
    #run_nf(aln, run_name)
    
    """
    Parse and save iteration outputs to the Results dict
    """
    bic, subtrees, alns = parse_outputs(aln, run, new_trees)
    temp_dict = {}
    temp_dict[run] = {
            "bic": bic,
            "input_trees": subtrees, 
            "aln": alns
            }
    print(temp_dict)
    return(temp_dict)

if __name__ == '__main__':
    """
    Initialise variables and Results dictionary
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    aln=args.aln
    tree_names=None
    Results=None
    n_iter=0

    """
    Take an input alignment, construct two subtrees from +R2 site assignment,
    input subtrees to MAST
    """
    Results_test = run_iteration(aln, tree_names, n_iter, Results)
