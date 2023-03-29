from Bio import SeqIO

class PartitioningScheme():
    """
    A PartitioningScheme instance retains summary information of the best 
    performing HMMSTER model retreived from a SubAlignment instance.

    Here, subaln == SubAlignment to easier distinguish from self.alignment
    """
    def __init__(self, alignment):
        self.alignment = alignment
        self.partitions = self.sites_from_alignment()
        self.bic_2t = float()
        self.model = ""

    def alignment_from_sites(self, subaln):
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

    def sites_from_alignment(self):
        """
        Generate a partition sites from an alignment
        i.e. reverse of alignment_from_sites

        As the goal here is to create a tuple (start, end) for each contiguous
        partition, the function will go through each site in the alignment and
        check if the partition at site i+1 differs from site i.
        """
        # Prepare dictionary of all possible partitions
        d = {key: [] for key in set(self.alignment)}
        
        current_partition = self.alignment[0]
        start_pos = 1
       
        for i, partition in enumerate(self.alignment):
            if partition != current_partition:
                # usually, the site position is the list index i+1 but we need
                # the position of the last partition before it changes,
                # therefore use i.
                prev_partition = self.alignment[i-1]
                end_pos = i
                # make tuple and append to dict value
                d[prev_partition].append((start_pos, end_pos))
                # Update to prepare for next partitions
                current_partition = partition
                start_pos = i+1
            if i+1 == len(self.alignment):
                # If at the end of the alignment, return end position and stop
                end_pos = i+1
                assert end_pos == len(self.alignment)
                prev_partition = self.alignment[i-1]
                d[prev_partition].append((start_pos, end_pos))
        return d
