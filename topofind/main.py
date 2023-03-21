import argparse
from topofind import partitioningscheme
from topofind import utils

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, 
            help='Sequence alignment (fasta)', required=True)
    return parser.parse_args()

def main():
    args = set_args()
    bic_improving=True

    """
    Initialise dictionaries for classes
    """
    # TODO: make dict to save all PartitioningSchemes
    # TODO: make dict to save all SubAlignments
    
    """
    Initialise first PartitioningScheme class
    """
    seq_len = utils.get_aln_length(args.aln)
    init_aln = utils.init_alignment(seq_len)
    PartitioningScheme = partitioningscheme.PartitioningScheme(init_aln, float("inf"))
    print(PartitioningScheme.alignment)
    print(PartitioningScheme.bic)
    while bic_improving:
        bic_improving=False

