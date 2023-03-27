import unittest
import os 
from topofind.partitioningscheme import PartitioningScheme
from topofind.subalignment import SubAlignment

class TestPartitioningScheme(unittest.TestCase):

    def setUp(self):
        """
        Create a dummy object for testing
        """
        self.subaln_test1 = SubAlignment()
        self.subaln_test1.sites_A = [(1,3500)]
        self.subaln_test1.sites_B = [(3501,5000)]
    
    def test_site_in_range(self):
        """
        Check that sites are assigned to the correct (discontiguous) partition,
        represented by tuples (start, end).
        """
        sites_A = [(1,50), (100,150)]
        sites_B = [(51,99), (151,200)]
        self.assertTrue(PartitioningScheme.site_in_range(self, 20, sites_A))
        self.assertFalse(PartitioningScheme.site_in_range(self, 20, sites_B))
        self.assertFalse(PartitioningScheme.site_in_range(self, 151, sites_A))
        self.assertTrue(PartitioningScheme.site_in_range(self, 151, sites_B))

    def test_create_alignment(self):
        self.partscheme_test1 = PartitioningScheme(self.subaln_test1)
        self.assertTrue(len(self.partscheme_test1.alignment) == 5000)
        self.assertTrue(self.partscheme_test1.alignment == ['A']*3500 + ['B']*1500)

if __name__ == '__main__':
    unittest.main()
