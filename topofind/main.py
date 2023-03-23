import argparse
import os
from topofind import partitioningscheme
from topofind.subalignment import SubAlignment
from topofind import utils

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln_path', type=str, 
            help='Sequence alignment (fasta)', required=True)
    parser.add_argument('-nt', '--num_threads', type=int, 
            help='Number of threads for IQ-TREE.', default=1, required=False)
    return parser.parse_args()

def main():
    repo_path=os.path.dirname(__file__)
    args=set_args()
    bic_improving=True
    n_trees=1

    """
    Initialise dictionaries for classes
    """
    # TODO: make dict to save all PartitioningSchemes
    # TODO: make dict to save all SubAlignments
    
    """
    Initialise first PartitioningScheme class
    """

    #seq_len = t1.get_aln_length(args.aln)
    #init_aln = t1.init_alignment(seq_len)
    #PartitioningScheme = partitioningscheme.PartitioningScheme(init_aln, float("inf"))
    
    while bic_improving:
        # 1. Get the different partitions from PartitioningScheme.alignment
        if n_trees == 1:
            # Create a SubAlignment instance and run R2 analysis
            # An instance is an individual object created from a class,
            # where a class is like a blueprint
            subalignment = SubAlignment()
            subalignment.run_r2(args.aln_path, args.num_threads)
            subalignment.run_Rhmm(repo_path)

        # 2. Identify the partition with the lowest BIC after previous split()
        # 3. Run split()
        # 4. Record results to SubAlignment and PartitioningScheme
        bic_improving=False

