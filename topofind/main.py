import argparse
import os
from topofind.partitioningscheme import PartitioningScheme
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
    Initialise dictionaries for result objects
    """
    partitioning_schemes = {}
    subalignments = {}

    while bic_improving:
        if not partitioning_schemes:
            # Create blank, starting alignment on first pass
            alignment = utils.new_alignment(args.aln_path)
            new_PartScheme = PartitioningScheme(alignment)

        possible_schemes = set(new_PartScheme.alignment)
        for partition_name in possible_schemes:
            SubAln = SubAlignment(partition_name)
            SubAln.iteration(args.aln_path, args.num_threads, repo_path) 
        bic_improving=False
