import unittest
from collections import OrderedDict
from run import compare_bic

class testCompareBIC(unittest.TestCase):

    def test_test1(self):
        '''
        Based on true negative test1.fa where t=2 BIC should be better than t=3
        '''
        test1=OrderedDict()
        test1['2_mast_A_B'] = {'bic': 81891.4428}
        test1['3_mast_B_AA_AB'] = None
        test1['3_mast_A_BA_BB'] = {'bic': 81972.3564}
        n_trees=3
        self.assertFalse(compare_bic(test1, n_trees))

    def test_no_runs(self):
        ''' what happens if no mast is run at the last iteration? '''
        no_runs=OrderedDict()
        no_runs['2_mast_A_B'] = {'bic': 81891.4428}
        no_runs['3_mast_B_AA_AB'] = None
        no_runs['3_mast_A_BA_BB'] = None
        n_trees=3
        self.assertFalse(compare_bic(no_runs, n_trees))

if __name__ == '__main__':
    unittest.main()
