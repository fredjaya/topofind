from Bio import SeqIO

class PartitioningScheme():
    """
    A PartitioningScheme instance retains summary information of the best 
    performing HMMSTER model retreived from a SubAlignment instance.

    Here, subaln == SubAlignment to easier distinguish from self.alignment
    """
    def __init__(self, subaln):
        self.alignment = self.create_alignment(subaln)
        self.bic_2t = subaln.bic_2t
        self.model = subaln.model

    def create_alignment(self, subaln):
        """
        Create a list where each element corresponds to a site.
        Each site corresponds to, and will be iteratively updated to
        the best-fit topology out of HMMSTER.

        This will be the core data/variable to access with subsequent splits.
        """
        alignment = []
        for i in range(1, subaln.total_sites()+1):
            if self.site_in_range(i, subaln.sites_A):
                alignment.append("A")
            elif self.site_in_range(i, subaln.sites_B):
                alignment.append("B")
        return alignment

    def site_in_range(self, site, ranges):
        """
        Determine which partition a site belongs to. 
        ranges a list of tuples (site_A or B)
        """
        for start, end in ranges:
            if start <= site <= end:
                return True
        return False

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
        return max(seq_lengthsd)

    def init_alignment(seq_length):
        return ["0"]*seq_length
