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
    parser.add_argument('-a', '--aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-f', '--aln_format', choices=['fasta', 'phylip'], help='Alignment format', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-T', '--num_threads', type=int, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    parser.add_argument('-r', '--run_nf', action='store_false', help='[DEBUG] Run nextflow processes?')
    parser.add_argument('-v', '--verbose', action='store_false', help='Print detailed output')
    args = parser.parse_args()
    return(args)

def run_nf(aln, run_name, mode, trees, submodel, n_trees):
    if args.run_nf:
        if n_trees == 2:
            aln_format=args.aln_format
        else:
            aln_format="fasta"
        nf = subprocess.run(["nextflow", "run", f"{repo_path}/python.nf",
            "--aln", aln,
            "--aln_format", aln_format,
            "--run_name", run_name,
            "--out", args.output_dir,
            "--nthreads", str(args.num_threads),
            "--mode", mode, 
            "--trees", trees, 
            "--submodel", submodel, 
            "-profile", args.executor],
            capture_output=args.verbose)
    else:
        return

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

def split_aln(aln, n_trees, tree_names, PartitionedTrees):
    run_name=f"{n_trees}_split_{'_'.join(tree_names)}"
    temp_out=f"{args.output_dir}/{run_name}"
    aln_name=os.path.basename(aln)

    print(f"[{run_name}]\tAssigning sites in {aln_name} to +R2 rate categories and making trees for each partition.")
    print(f"[{run_name}]\tOutput trees: {tree_names}")

    run_nf(aln, run_name, "split_aln", "null", "null", n_trees) 

    '''Save trees to dict'''
    trees=sorted(glob.glob(f"{temp_out}/*-out.treefile"))
    for key, val in zip(tree_names, trees):
        PartitionedTrees[key] = val
    return

    print(f"[{run_name}]\tDone! Files output to {temp_out}")

def n_models(n_trees):
    return(",".join(["GTR+FO+G"]*n_trees))

def mast(n_trees, tree_names, PartitionedTrees):
    run_name=f"{n_trees}_mast_{'_'.join(tree_names)}"
    temp_out=f"{args.output_dir}/{run_name}"
    concat_tree=f"{temp_out}/concat.treefile"
    if not os.path.exists(temp_out):
        os.mkdir(temp_out)
    
    ''' Concatenate n trees to input to MAST '''
    to_concat = []
    for key, value in PartitionedTrees.items():
        if key in tree_names:
           to_concat.append(value)

    ''' Check that there is a new tree for MAST, else don't run. Can refactor '''
    if None in to_concat:
        print(f"\033[1m[{run_name}]\tNo new partitions, MAST not run.\033[0m")
        MastResults[run_name] = None

    elif to_concat == []:
        """ 
        test2.fa: ERROR: Only one state is observed in alignment
        Error handled by check_valid_runs()
        """
        print(f"\033[1m[{run_name}]\tNo new partitions, MAST not run.\033[0m")
        MastResults[run_name] = None

    elif len(to_concat) == 1:
        """
        test3.fa: ERROR: All sites assigned to a single tree
        """
        print(f"\033[1m[{run_name}]\tNo new partitions, MAST not run.\033[0m")
        MastResults[run_name] = None

    else:
        with open(concat_tree, "w") as file:
            subprocess.run(["cat"] + to_concat, stdout=file)

        print(f"[{run_name}]\tRunning MAST with Trees: {tree_names} as input.")
        run_nf(args.aln, run_name, "mast", concat_tree, n_models(n_trees), n_trees)
        print(f"[{run_name}]\tDone! Files output to {temp_out}")

        """
        Parse and save iteration outputs to the Results dict
        """
        
        '''Record MAST BIC'''
        iqtree_out_path=f"{temp_out}/t2_mast_tr.iqtree"
        bic=get_bic(iqtree_out_path)
        print(f"\033[1m[{run_name}]\tBIC: {bic}\033[0m")

        '''Collect fastas post-HMM assignment'''
        alns=collect_alns(glob.glob(f"{temp_out}/*.fas"), tree_names)
   
        MastResults[run_name] = {
            "bic": bic,
            "input_trees": tree_names, 
            "aln": alns
            }
        return

def compare_bic(MastResults, n_trees):
    temp_dict={}
    for key, value in MastResults.items():
        try:
            temp_dict[key.split('_')[0]] = value["bic"]
        except TypeError:
            """ When mast isn't run, value["bic"] == None """
            pass
    print("\nHas the BIC improved with more trees?")
    best_n=min(temp_dict, key=temp_dict.get)
    if int(best_n) < n_trees:
        print(f"\nNo improvement in BIC. Stopping program :)")
        return False
    else:
        print(f"\nYes! Proceed with more trees.")
        return True

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
    Initialise variables and Results dictionary
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    n_trees=2
    tree_names=['A', 'B']
    bic_improving=True
    PartitionedTrees=OrderedDict()
    MastResults=OrderedDict()
    
    print(f"\n\
        aln:            {args.aln}\n\
        aln_format:     {args.aln_format}\n\
        output_dir:     {args.output_dir}\n\
        num_threads:    {args.num_threads}\n\
        nf_executor:    {args.executor}\n"
        )

    """
    Always run first iteration (2 trees --> MAST)
    """
    split_aln(args.aln, n_trees, tree_names, PartitionedTrees)
    mast(n_trees, tree_names, PartitionedTrees)

    """
    Always run second iteration (3 trees --> MAST)
    """
    while bic_improving:
        check_valid_runs(MastResults)
        n_trees+=1
        for t_old in tree_names: 
            # TODO: Run in parallel
            ''' Split each existing alignment and add key to existing dict'''
            p=[f"{t_old}A", f"{t_old}B"]
            for t_new in p:
                PartitionedTrees[t_new] = None
            ''' Access alignment by matching key from previous MAST run '''
            pattern = '_'.join(tree_names)
            aln=[value['aln'][f"{t_old}′"] for key, value in MastResults.items() if key.endswith(pattern)][0]

            split_aln(aln, n_trees, p, PartitionedTrees) 
        
        '''Get list of tree combinations for next iteration'''
        for i in range(0, len(tree_names)):
            '''Name the run'''
            tree_list=recurse_trees(tree_names, i)
            run_name=f"{n_trees}_mast_{'_'.join(tree_list)}"
            '''Add required input trees'''
            MastResults[run_name] = {"input_trees": tree_list}
            mast(n_trees, tree_list, PartitionedTrees)
            bic_improving = compare_bic(MastResults, n_trees)

    """
    Exiting program
    """
    print(f"\nFinal Results:")
    print("\n", PartitionedTrees, "\n")
    print(MastResults)
    
    """
    TODO:
        Do not run split_aln if already exists. 
        Alternatively, remake the tree using output MAST partitions.
        For example, For [A, BA, BB] --> MAST, Tree [A] and Alignment [A] already exists. 
        However, Alignment [A`] exists, but no tree is constructed from [A`].
    """

