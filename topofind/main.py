import argparse

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln', type=str, help='Sequence alignment (fasta)', required=True)
    return parser.parse_args()

def main():
    args = set_args()

    return
