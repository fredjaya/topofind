import uuid
import os
import re 
from topofind import utils
from Bio import SeqIO
from io import StringIO

class SubAlignment():
    def __init__(self):
        self.rid = self.random_id()
        self.model = ""
        #self.r2_partition_A = {} # k: header; v: seq
        #self.r2_partition_B = {} 
        self.bic_1t = float()
        self.bic_2t = float()
        self.topology_A = ""
        self.topology_B = ""
        self.sites_A = []
        self.sites_B = []

    def total_sites(self):
        return max(self.sites_A[-1][1], self.sites_B[-1][1])

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
        print(f"[{self.rid}]\tMaking a single tree with the best-fitting +R2 model")
        os.mkdir(self.rid)
        # TODO: Replace hardcode
        iqtree_path="/home/frederickjaya/Downloads/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
        cmd = f"{iqtree_path} -s {aln_path} -pre {self.rid}/r2 -mrate R2 -nt {str(nthreads)} -wslr -wspr -alninfo"
        stdout, stderr, exit_code = utils.run_command(cmd)

    def parse_best_model(self):
        """
        "grep" the best substitution model from run_r2 and save as attribute
        for downstream iq-tree runs
        """
        iqtree_file = f"{self.rid}/r2.iqtree"
        r = re.compile("^Best-fit model according to BIC: ")

        with open(iqtree_file, 'r') as f:
            for line in f:
                if re.search(r, line):
                    self.model = line.split(" ")[-1].strip()
                    print(f"[{self.rid}]\tBest-fit model according to BIC: {self.model}")
    
    def parse_bic(self, iqtree_file, attr):
        """
        "grep" the BIC value from an .iqtree file and 
        save to either bic_1t or bic_2t
        """
        r = re.compile("\(BIC\)")
        with open(iqtree_file, 'r') as f:
            for line in f:
                if re.search(r, line):
                    line = float(line.split(" ")[-1].strip())
                    setattr(self, attr, line)

    def run_Rhmm(self, repo_path):
        """
        Using the MixtureModelHMM in R, assign sites to one of two +R2 classes
        """
        print(f"[{self.rid}]\tAssigning sites to rate-categories with HMM")
        rscript_path = os.path.join(repo_path, "../bin", "hmm_assign_sites.R")
        sitelh = f"{self.rid}/r2.sitelh"
        alninfo = f"{self.rid}/r2.alninfo"
        cmd = f"Rscript {rscript_path} {sitelh} {alninfo} {self.rid}"
        stdout, stderr, exit_code = utils.run_command(cmd)

    def partition_aln(self, aln_file):
        """
        Partition the previous alignment into separate .fa files  

        TODO: 
            - check if previous alignment exists in memory
        """
        print(f"[{self.rid}]\tPartitioning...")
        # Read partitions and save to list
        partitions = []
        part_path = f"{self.rid}/r2.partition"
        with open(part_path, 'r') as pfile:
            for line in pfile:
                if line.startswith("\tcharset"):
                    # Convert to tuple and save
                    start, end = re.findall(r"(\d+-\d+)", line)[0].split("-")
                    start, end = int(start), int(end)
                    partitions.append((start, end))

        # Read alignment to split and store as dict
        sequences = {}
        for seq in SeqIO.parse(aln_file, "fasta"):
            sequences[seq.id] = seq
        
        # Write partitioned fasta files
        # enumerate iterates over an objects index and element
        for i, partition in enumerate(partitions):
            # Each alignment is always split according to R2
            pname = "A"
            if i == 1:
                pname = "B"
            start, end = partition
            out_name = f"{self.rid}/partition_{pname}.fasta"

            with open(out_name, "w") as ofile:
                for header, whole_seq in sequences.items():
                    # start-1 because 0-based indexing
                    part_seq = str(whole_seq.seq[start-1:end])
                    ofile.write(f">{header}\n")
                    ofile.write(f"{part_seq}\n")

            """
            for header, seq_record in sequences.items():
                partitioned_seq = str(seq_record.seq[start-1:end])
                if pname == 'A':
                    self.r2_partition_A[header] = partitioned_seq
                elif pname == 'B':
                    self.r2_partition_B[header] = partitioned_seq
            """

    def run_iqtree_on_parts(self, pname, nthreads):
        """
        Run iqtree using the best model from run_r2 on each new partition.
        """
        print(f"[{self.rid}]\tMaking trees for partition {pname}")
        file_prefix = f"{self.rid}/partition_{pname}"
        # TODO: Replace hardcode
        iqtree_path="/home/frederickjaya/Downloads/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
        cmd = f"{iqtree_path} -s {file_prefix}.fasta -pre {file_prefix} -m {self.model} -nt {nthreads}"
        stdout, stderr, exit_code = utils.run_command(cmd)

    def concat_part_trees(self):
        """
        Concatenate trees made from part A and B into a single file for tree mixture input
        """
        # TODO: Account for only one topology
        with open(f"{self.rid}/partition_A.treefile", 'r') as tfile_A, \
        open(f"{self.rid}/partition_B.treefile", 'r') as tfile_B, \
        open(f"{self.rid}/concat.treefile", 'w') as outfile:
            outfile.write(tfile_A.read())
            outfile.write(tfile_B.read())

    def run_hmmster(self, original_aln, num_threads):
        print(f"[{self.rid}]\tRunning HMMSTER")
        iqtree_path="/home/frederickjaya/Downloads/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
        cmd = f"{iqtree_path} -s {original_aln} -pre {self.rid}/hmmster -m {self.model}+T -nt {num_threads} -hmmster{{gm}} -te {self.rid}/concat.treefile"
        stdout, stderr, exit_code = utils.run_command(cmd)

    def parse_topologies(self):
        """
        Store output HMMSTER topologies
        """
        with open(f"{self.rid}/hmmster.treefile", 'r') as f:
            for line, pname in zip(f, ['A', 'B']):
                setattr(self, f"topology_{pname}", line.strip())

    def parse_hmm_sites(self):
        """
        Store output HMMSTER topologies
        """
        r = re.compile(r"\[\d+,\d+\]\t\d+$")
        with open(f"{self.rid}/hmmster.hmm", 'r') as f:
            for line in f:
                if line.startswith("["):
                    line = line.split("\t")
                    start = int(re.findall("(?<=\[)\d+(?=,)", line[0])[0])
                    end = int(re.findall("(?<=,)\d+(?=\])", line[0])[0])
                    category = line[1].strip()

                    if category == "1":
                        self.sites_A.append((start, end))
                    elif category == "2":
                        self.sites_B.append((start, end))

    def iteration(self, aln, num_threads, repo_path):
        """
        Main pipeline for a single iteration of +R2, splitting, and HMMSTER
        """
        # TODO: What if R1 > R2?
        self.run_r2(aln, num_threads)
        self.parse_best_model()
        self.parse_bic(f"{self.rid}/r2.iqtree", "bic_1t")
        self.run_Rhmm(repo_path)
        self.partition_aln(aln)
        # TODO: Run in parallel
        self.run_iqtree_on_parts("A", num_threads)
        self.run_iqtree_on_parts("B", num_threads)
        self.concat_part_trees()
        self.run_hmmster(aln, num_threads)
        # TODO: Get correct HMM BIC
        self.parse_bic(f"{self.rid}/hmmster.iqtree", "bic_2t")
        self.parse_topologies() 
        self.parse_hmm_sites() 
        print(vars(self))
