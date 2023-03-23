import uuid
import os
import subprocess 

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
        iqtree_path="/home/frederickjaya/Downloads/iqtree-2.2.3.hmmster-Linux/bin/iqtree2"
        cmd = f"{iqtree_path} -s {aln_path} -pre {self.rid}/r2 -mrate R2 -nt {str(nthreads)} -wslr -wspr -alninfo"
        print(cmd)
        # Popen lets you access the I/O pipes. stdout and stderr specify which O pipes you want to access.
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        exit_code = process.returncode
    
        if exit_code != 0:
            print("Error")
