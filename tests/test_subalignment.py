import unittest
import os 
from topofind.subalignment import SubAlignment

class TestSubAlignment(unittest.TestCase):
    # This can be placed in the class-level (instead of in setUp)
    # because it does not need  to be reinitialised before every test
    repo_path=os.path.dirname(__file__)

    def test_iteration_test1(self):
        subaln = SubAlignment()
        aln = f"{self.repo_path}/../data/test1.fa"
        subaln.iteration(aln, 4, self.repo_path)

        self.assertTrue(subaln.model == 'JC+R2')
        self.assertAlmostEqual(subaln.bic_1t, 83581.5, delta=0.1)
        self.assertAlmostEqual(subaln.bic_2t, 81088.3, delta=0.1)
        self.assertTrue(subaln.sites_A == [(1,3500)])
        self.assertTrue(subaln.sites_B == [(3501,5000)])

if __name__ == '__main__':
    unittest.main()
