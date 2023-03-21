from Bio import SeqIO

def get_aln_length(aln_path):
    """
    Get the longest sequence length to initialise the first 
    PartitioningScheme.alignment = []
    """
    alns = SeqIO.parse(aln_path, "fasta")
    seq_lengths = set()
    for a in alns:
        seq_lengths.update({len(a.seq)})
    # TODO: What if unaligned, or length is variable?
    return max(seq_lengths)

def init_alignment(seq_length):
    return ["0"]*seq_length
