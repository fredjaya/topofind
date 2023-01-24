import argparse
import os
import subprocess
import re
import glob
import multiprocessing 
from collections import OrderedDict

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

def pool_iterations(aln, tree_names, n_iter, Results):
    with multiprocessing.Pool() as pool:
        pool.map(run_iteration, aln, tree_names, n_iter, Results)

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

def recurse_trees(tree_list, pos):
    """
    This operates on a list of existing trees and makes a new
    list of trees for the next iteration
    """
    if tree_list == None:
        return(['A', 'B'])
    else:
        try:
            new_list=tree_list.copy()
            t0 = str(new_list[pos])
            t1 = f"{t0}A" 
            t2 = f"{t0}B"
            new_list += [t1, t2]
            new_list.remove(t0)
            return(new_list)
        except AttributeError:
            '''When tree_list==None, you cant .copy()'''
            pass

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
    new_trees = recurse_trees(tree_names, 0)
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
    Results[run] = {
            "bic": bic,
            "input_trees": subtrees, 
            "aln": alns
            }
    return(Results, run, n_iter, new_trees)

def get_n_last_runs():
    """
    Calculate how many combinations of trees were input to MAST
    in the last iteration
    """

if __name__ == '__main__':
    """
    Initialise variables and Results dictionary
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    aln=args.aln
    n_iter=0
    tree_names=None
    Results=OrderedDict()

    """
    Take an input alignment, construct two subtrees from +R2 site assignment,
    input subtrees to MAST
    """
    if n_iter < 2 : 
        '''Always run two iterations'''
        Results, prev_runs, n_iter, tree_names = run_iteration(aln, tree_names, n_iter, Results)
        #n_last_runs=
        #last_runs=list(Results.keys())[n_last_runs]
        '''Get list of tree combinations for next iteration'''
        n_iter+=1
        for i in range(0, len(tree_names)):
            '''Name the run'''
            tree_list=recurse_trees(tree_names, i)
            run=f"{n_iter}_{'_'.join(tree_list)}"
            Results[run] = {}
            '''Add required input trees'''
            Results[run]["input_trees"] = tree_list
        print(Results)
        '''Get output alignments/partitions from previous iteration'''
        pool_cmds=[
            list(Results[prev_runs]['aln'].values()), 
            list(Results[prev_runs]['aln'].keys()),
            n_iter, 
            Results
            ]

        #print(pool_cmds)
        '''with multiprocessing.Pool() as pool:
            pool.map(run_iteration, pool_cmds)'''
        #run_iteration(list(Results[prev_runs]['aln'].values())[0], list(Results[prev_runs]['aln'].keys())[0], n_iter, Results)

    else:
        print("else")
