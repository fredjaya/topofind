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

    def test_sites_from_alignment(self):
        """
        Simple two-partition case
        """
        # Two-partition case
        PS1 = PartitioningScheme(['A', 'A', 'A', 'B', 'B', 'A', 'A'])
        self.assertEqual(PS1.partitions, {'A': [(1,3),(6,7)], 'B': [(4,5)]})
    
        """
        Three partitions
        """
        PS2 = PartitioningScheme(['AA', 'AA', 'AB', 'B', 'B', 'AB', 'AB'])
        self.assertEqual(PS2.partitions, {'AA': [(1,2)], 'AB': [(3,3),(6,7)], 'B': [(4,5)]})

        """
        Three partitions with single-site partitions
        """
        PS3 = PartitioningScheme(['AA', 'AA', 'ABA', 'B', 'B', 'ABA', 'ABB'])
        self.assertEqual(PS3.partitions, {'AA': [(1,2)], 'ABA': [(3,3),(6,6)], 'ABB': [(7,7)], 'B': [(4,5)]})
if __name__ == '__main__':
    unittest.main()
