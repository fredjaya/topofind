import uuid
import os
import re 
from topofind import utils
from Bio import SeqIO

class SubAlignment():
    def __init__(self):
        self.rid = self.random_id()
        self.bic_1t = float()
        self.bic_2t = float()
        self.topo_A = ""
        self.topo_B = ""
        self.sites_A = []
        self.sites_B = []

    def sites_all(self):
        # Moved out of __init__ as it can be derived from existing attributes
        return self.sites_A + self.sites_B

    def random_id(self):
        """
        Generate a unique (hopefully) string to identify each run
        """
        return str(uuid.uuid4()).split("-")[0]
    
    def run_r2(self, aln_path, nthreads):
        """
        Take an input alignment, construct a single tree using the best-fitting +R2
        model and output site-lh and alninfo files.
        """
        os.mkdir(self.rid)
        # TODO: Replace hardcode
        iqtree_path="/home/frederickjaya/Downloads/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
        cmd = f"{iqtree_path} -s {aln_path} -pre {self.rid}/r2 -mrate R2 -nt {str(nthreads)} -wslr -wspr -alninfo"
        stdout, stderr, exit_code = utils.run_command(self.rid, cmd)
        if exit_code != 0:
            # TODO: Deal with different cases.
            print(stdout)

    def run_Rhmm(self, repo_path):
        """
        Using the MixtureModelHMM in R, assign sites to one of two +R2 classes
        """
        rscript_path = os.path.join(repo_path, "../bin", "hmm_assign_sites.R")
        sitelh = f"{self.rid}/r2.sitelh"
        alninfo = f"{self.rid}/r2.alninfo"
        cmd = f"Rscript {rscript_path} {sitelh} {alninfo} {self.rid}"
        stdout, stderr, exit_code = utils.run_command(self.rid, cmd)
        if exit_code != 0:
            # TODO: Deal with different cases.
            print(stdout)

    def partition_aln(self, aln_file):
        """
        Partition the previous alignment into separate .fa files  

        TODO: 
            - check if previous alignment exists in memory
        """
        # Read partitions and save to list
        partitions = []
        part_path = f"{self.rid}/r2.partition"
        with open(part_path, 'r') as pfile:
            for line in pfile:
                if line.startswith("\tcharset"):
                    partitions.append(re.findall(r"(\d+-\d+)", line)[0])
        print(partitions)

        # Read sequences and store as dict
        sequences = {}
        for seq in SeqIO.parse(aln_file, "fasta"):
            sequences[seq.id] = seq
        
        # Write partitioned fasta files
        # enumerate iterates over an objects index and element
        for i, partition in enumerate(partitions):
            pname = "A"
            if i == 1:
                pname = "B"
            start, end = map(int, partition.split("-"))
            out_name = f"{self.rid}/partition_{pname}.fasta"

            # TODO: what if partitions are not contiguous
            with open(out_name, "w") as ofile:
                for header, seq_record in sequences.items():
                    ofile.write(f">{header}\n")
                    ofile.write(f"{str(seq_record.seq[start-1:end])}\n")
