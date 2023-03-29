import subprocess
from Bio import SeqIO

def run_command(cmd):
    """
    General function to run command-line processes.

    Usage:
        stdout, stderr, exit_code = utils.run_command(cmd)
    """
    # Popen lets you access the I/O pipes.
    # stdout and stderr options specify which pipes you want to capture.
    process = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # communicate() waits for the process to complete and returns a tuple of
    # the captured pipes.
    stdout, stderr = process.communicate()
    exit_code = process.returncode
    # Print error message
    if exit_code != 0:
        print(stderr)
    return stdout, stderr, exit_code

def new_alignment(aln_path):
    """                                                    
    Get the longest sequence length to initialise the first
    PartitioningScheme.alignment = []                      
    """
    alns = SeqIO.parse(aln_path, "fasta")
    seq_lengths = set()                  
    for a in alns:                       
        seq_lengths.update({len(a.seq)}) 
    return ["0"]*max(seq_lengths)        
