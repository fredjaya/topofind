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

def run_nf(aln, run_name, mode, trees, submodel):
    subprocess.run(["nextflow", "run", f"{repo_path}/python.nf",
        "--aln", aln,
        "--aln_format", args.aln_format,
        "--run_name", run_name,
        "--out", args.output_dir,
        "--nthreads", str(args.num_threads),
        "--mode", mode, 
        "--trees", trees, 
        "--submodel", submodel, 
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

    '''Record MAST BIC'''
    iqtree_out_path=f"{out_path}/t2_mast_tr.iqtree"
    bic=get_bic(iqtree_out_path)

    '''Collect subtrees used as input to mast'''
    subtrees=collect_subtrees(f"{out_path}/concat.treefile", new_trees)

    '''Collect fastas post-HMM assignment'''
    ## Get partitions
    ## Get alignment files
    alns=collect_alns(glob.glob(f"{out_path}/*.fas"), new_trees)
    return(bic, subtrees, alns)

def split_aln(aln, n_iter, tree_names, PartitionedTrees):
    run_name=f"{n_iter}_split_{'_'.join(tree_names)}"
    temp_out=f"{args.output_dir}/{run_name}"
    aln_name=os.path.basename(aln)

    print(f"[{run_name}]\tAssigning sites in {aln_name} to +R2 rate categories and making trees for each partition.")
    print(f"[{run_name}]\tOutput trees: {tree_names}")

    run_nf(aln, run_name, "split_aln", "null", "null") 
    print(f"[{run_name}]\tDone! Files output to {temp_out}")

    '''Save trees to dict'''
    trees=sorted(glob.glob(f"{temp_out}/*-out.treefile"))
    for key, val in zip(tree_names, trees):
        PartitionedTrees[key] = val
    return

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
    n_iter=1
    tree_names=['A', 'B']
    PartitionedTrees=OrderedDict()
    MastResults=OrderedDict()

    print(f"\n\
        aln:            {args.aln}\n\
        aln_format:     {args.aln_format}\n\
        output_dir:     {args.output_dir}\n\
        num_threads:    {args.num_threads}\n\
        nf_executor:    {args.executor}\n"
        )

    split_aln(args.aln, n_iter, tree_names, PartitionedTrees)
    '''mast time'''
    '''First, concatenate tree files'''

    run_name=f"{n_iter}_mast_{'_'.join(tree_names)}"
    temp_out=f"{args.output_dir}/{run_name}"
    concat_tree=f"{temp_out}/concat.treefile"
    if not os.path.exists(temp_out):
        os.mkdir(temp_out)
    
    with open(concat_tree, "w") as file:
        subprocess.run(["cat"] + list(PartitionedTrees.values()), stdout=file)

    print(f"[{run_name}]\tRunning MAST with Trees: {tree_names} as input.")
    #run_nf(args.aln, run_name, "mast", concat_tree, "GTR+FO+G,GTR+FO+G")
    print(f"[{run_name}]\tDone! Files output to {temp_out}")

    """
    Parse and save iteration outputs to the Results dict
    """
    
    '''Record MAST BIC'''
    iqtree_out_path=f"{temp_out}/t2_mast_tr.iqtree"
    bic=get_bic(iqtree_out_path)
    print(f"[{run_name}]\tBIC: {bic}")

    '''Collect fastas post-HMM assignment'''
    ## Get partitions
    ## Get alignment files
    alns=collect_alns(glob.glob(f"{temp_out}/*.fas"), tree_names)
   
    MastResults[run_name] = {
        "bic": bic,
        "input_trees": tree_names, 
        "aln": alns
        }
    print(f"[{run_name}]\t{PartitionedTrees}")
    print(f"[{run_name}]\t{MastResults}")
    """
    Run split_aln on one alignment
    """

    n_iter+=1
    for t_old in tree_names: 
        ''' Split each existing alignment and add key to existing dict'''
        p=[f"{t_old}A", f"{t_old}B"]
        run_name=f"{n_iter}_split_{'_'.join(p)}"
        for t_new in p:
            PartitionedTrees[t_new] = None
        ''' Access alignment by matching key from previous MAST run '''
        pattern = '_'.join(tree_names)
        aln=[value['aln'][f"{t_old}′"] for key, value in MastResults.items() if key.endswith(pattern)][0]
        print(aln[0])
        split_aln(aln, n_iter, p, PartitionedTrees)
    print(PartitionedTrees)
    """
    return(Results, run_name, n_iter, new_trees)
    Results, prev_runs, n_iter, tree_names = run_first_iter(aln, tree_names, n_iter, Results)
    print(Results)
    print(prev_runs)
    '''Run split_aln on one output alignment'''
    print(list(Results[prev_runs]['aln'])[0])

    for aln_name, file in Results[prev_runs]['aln'].items():
        print(f"\nSplitting partition {aln_name}")
        print(file)
        run_nf(
            aln = file,
            run_name = f"{prev_runs}/{aln_name}",
            mode = "split_aln",
            trees = "null",
            submodel = "null"
        )
    '''Run split_aln on both output alignments in parallel'''

    '''Get list of tree combinations for next iteration'''
    n_iter+=1
    for i in range(0, len(tree_names)):
        '''Name the run'''
        tree_list=recurse_trees(tree_names, i)
        run=f"{n_iter}_{'_'.join(tree_list)}"
        Results[run] = {}
        '''Add required input trees'''
        Results[run]["input_trees"] = tree_list
    """
