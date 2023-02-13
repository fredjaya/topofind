import argparse
import os
import subprocess
import re
import glob
import multiprocessing as mp
from collections import OrderedDict
import sys

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, help='Input sequence alignment', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Path to output directory', default=os.getcwd())
    parser.add_argument('-T', '--num_threads', type=int, help='[IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)', default=1)
    parser.add_argument('-e', '--executor', choices=['local', 'slurm'], help='[NEXTFLOW] Where to run nextflow processes (default: local)', default='local')
    parser.add_argument('-v', '--verbose', action='store_false', help='Print detailed output')
    parser.add_argument('-r', '--run_nf', action='store_false', help='[DEBUG] Run nextflow processes?')
    parser.add_argument('-s', '--skip_iter', type=int, help='[DEBUG] Only run nextflow processes from n_trees', default=0)
    args = parser.parse_args()
    return(args)

def run_nf(aln, run_name, mode, trees, submodel, n_trees):
    if args.run_nf:
        if args.skip_iter < n_trees:
            nf = subprocess.run(["nextflow", "run", f"{repo_path}/python.nf",
                "--aln", aln,
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

def n_models(n_trees):
    return(",".join(["GTR+FO+G"]*n_trees))

def mast(n_trees, run_name, tree_names, PartitionedTrees):
    temp_out=f"{args.output_dir}/{run_name}"
    temp_dict = {}
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
        temp_dict[run_name] = None

    elif to_concat == []:
        """ 
        test2.fa: ERROR: Only one state is observed in alignment
        Error handled by check_valid_runs()
        """
        print(f"\033[1m[{run_name}]\tNo new partitions, MAST not run.\033[0m")
        temp_dict[run_name] = None

    elif len(to_concat) == 1:
        """
        test3.fa: ERROR: All sites assigned to a single tree
        """
        print(f"\033[1m[{run_name}]\tNo new partitions, MAST not run.\033[0m")
        temp_dict[run_name] = None

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
        
        temp_dict[run_name] = {
            "bic": bic,
            "input_trees": tree_names, 
            "aln": alns
            }
    return(temp_dict)

def get_new_trees(MastResults, n_trees):
    ''' Get previous input trees from last mast iterations '''
    tree_names = []
    for key, value in MastResults.items():
        if key.startswith(str(n_trees-1)):
            if value is not None:
                tree_names.append(value['input_trees'])
    ''' For each combination of input trees, make new ones by splitting one tree '''
    new_trees = []
    for l in tree_names:
        for i in range(0, len(l)):
            new_trees.append(recurse_trees(l, i))
    return(new_trees)
    
def compare_bic(MastResults, n_trees):
    temp=[]
    for key, value in MastResults.items():
        try:
            temp.append((key, value["bic"]))
        except TypeError:
            """ When mast isn't run, value["bic"] == None """
            pass
    print("\nHas the BIC improved with more trees?")
    best_model=min(temp, key = lambda x: x[1])
    if int(best_model[0].split("_")[0]) < n_trees:
        print(f"No improvement in BIC. Stopping program :)\n")
        return False
    else:
        print(f"Yes! Proceed with more trees.\n")
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

def new_partitions_per_tree(tree):
    """ 
    Each iteration requires all trees in a set be split.
    Input is a single tree accessed by looping through the list of trees
    e.g. 'A' from ['A', 'B'] --> [AA, AB]
    """
    return([f"{tree}A", f"{tree}B"])

def get_previous_run(n_trees, tree_list):
    """
    Return the key to access the previous alignment required for further partitioning
    e.g. to get new partitions for tree 'AA' or 'AB', alignment 'A`' from MastResults
    is needed. Alignment 'A`' exists in run 2_mast_A_B.
    """
    t='_'.join(tree_list)
    return(str(f"{n_trees-1}_mast_{t}"))

def get_previous_aln(prev_run, tree, MastResults):
    """
    With the run name (key) from get_previous_run(),
    now access MastResults for the aln path
    """
    for key, value in MastResults.items():
        if value is None:
            return None
        elif key in prev_run:
            return(value['aln'][f"{tree}′"])

if __name__ == '__main__':
    """
    Initialise variables and Results dictionary
    """
    args=set_args()
    repo_path=os.path.dirname(__file__)
    n_trees=2
    all_tree_names=[['A', 'B']]
    bic_improving=True
    PartitionedTrees=OrderedDict()
    MastResults=OrderedDict()
    
    print(f"\n\
        aln:            {args.aln}\n\
        output_dir:     {args.output_dir}\n\
        num_threads:    {args.num_threads}\n\
        nf_executor:    {args.executor}\n"
        )

    while bic_improving:
        if n_trees == 2:
            """ Always run first iteration """

            """ Set variables for split_aln and mast """
            input_trees = all_tree_names[0]
            run_name=f"{n_trees}_split_{'_'.join(input_trees)}"
            temp_out=f"{args.output_dir}/{run_name}"
            aln_name=os.path.basename(args.aln)
            print(f"[{run_name}]\tAssigning sites in {aln_name} to +R2 rate categories and making trees for each partition.")

            """ Run split_aln """
            run_nf(args.aln, run_name, "split_aln", "null", "null", n_trees) 

            """Save output trees to dict"""
            trees=sorted(glob.glob(f"{temp_out}/*-out.treefile"))
            for key, val in zip(input_trees, trees):
                PartitionedTrees[key] = val
            print(f"[{run_name}]\tDone! Files output to {temp_out}")

            """ Run mast """
            run_name=f"{n_trees}_mast_{'_'.join(input_trees)}"
            temp_dict = mast(n_trees, run_name, input_trees, PartitionedTrees)
            MastResults.update(temp_dict)

        check_valid_runs(MastResults)
        n_trees+=1
        print(f"\nRUNNING PROGRAM WITH {n_trees} TREES\n")
        
        for tree_list in all_tree_names:
            """ 
            Run split_aln to further partition all partitions from the previous iteration
            all_tree_names = [['A', 'BA', 'BB'], ['B', 'BA', 'BB']]
            tree_list = ['A', 'BA', 'BB']
            tree = 'A'
            """
            prev_run = get_previous_run(n_trees, tree_list)
            for tree in tree_list:
                """ 
                Further partition each existing alignment
                new_trees = ['AA', 'AB'] 
                new_tree  = 'AA'
                """
                new_trees=new_partitions_per_tree(tree)
                run_name=f"{n_trees}_split_{'_'.join(new_trees)}"
                temp_out=f"{args.output_dir}/{run_name}"
                aln=get_previous_aln(prev_run, tree, MastResults)
                
                run_split = True
                for new_tree in new_trees:
                    """ Run split_aln if it doesn't exist yet """
                    if new_tree in PartitionedTrees.keys():
                        run_split = False
                
                if run_split and aln is not None:
                    ''' Access previous alignment by matching key from previous MAST run '''
                    aln_name=os.path.basename(aln)
                    print(f"[{run_name}]\tAssigning sites in {aln_name} to +R2 rate categories and making trees for each partition.")
                    print(f"[{run_name}]\tOutput trees: {new_trees}")
                    run_nf(aln, run_name, "split_aln", "null", "null", n_trees) 
                    print(f"[{run_name}]\tDone! Files output to {temp_out}")

                '''Save trees to dict'''
                tree_files=sorted(glob.glob(f"{temp_out}/*-out.treefile"))
                if len(tree_files) != 2:
                    for new_tree in new_trees:
                        PartitionedTrees[new_tree]=None
                else:
                    for key, val in zip(new_trees, tree_files):
                        PartitionedTrees[key] = val
       
        print(PartitionedTrees)
        '''Get list of new possible tree combinations for MAST'''
        all_tree_names = get_new_trees(MastResults, n_trees)
        cmds=[]
        for tree_list in all_tree_names:
            ''' Run MAST on each combination of trees '''
            run_name=f"{n_trees}_mast_{'_'.join(tree_list)}"
            t = (n_trees, run_name, tree_list, PartitionedTrees)
            cmds.append(t)
        with mp.Pool() as pool:
            temp_dicts = pool.starmap(mast, cmds)
            for i in temp_dicts:
                MastResults.update(i)
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
