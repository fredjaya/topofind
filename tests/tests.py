import unittest
from collections import OrderedDict
from run import compare_bic, get_new_trees
from run import get_new_trees
from run import new_partitions_per_tree
from run import get_previous_run
from run import get_previous_aln

test1=OrderedDict()
test1['2_mast_A_B'] = {'bic': 81891.4428, 'input_trees': ['A', 'B'], 
        'aln': {'A′': 'test1_class_1-out.fas', 'B′': 'test1_class2-out.fas'}}
test1['3_mast_B_AA_AB'] = None
test1['3_mast_A_BA_BB'] = {'bic': 81972.3564, 'input_trees': ['A', 'BA', 'BB']}

no_runs=OrderedDict()
no_runs['2_mast_A_B'] = {'bic': 81891.4428}
no_runs['3_mast_B_AA_AB'] = None
no_runs['3_mast_A_BA_BB'] = None

imp=OrderedDict()
imp['2_mast_A_B'] = {'bic': 81891.4428, 'input_trees': ['A', 'B']}
imp['3_mast_B_AA_AB'] = {'bic': 80000, 'input_trees': ['B', 'AA', 'AB']}
imp['3_mast_A_BA_BB'] = {'bic': None}
imp['3_mast_A_BA_BB'] = {'bic': 75000, 'input_trees': ['A', 'BA', 'BB']}

class testCompareBIC(unittest.TestCase):

    def test_test1(self):
        '''
        Based on true negative test1.fa where t=2 BIC should be better than t=3
        '''
        n_trees=3
        self.assertFalse(compare_bic(test1, n_trees))

    def test_no_runs(self):
        ''' what happens if no mast is run at the last iteration? '''
        n_trees=3
        self.assertFalse(compare_bic(no_runs, n_trees))

    def test_improvement(self):
        n_trees=3
        self.assertTrue(compare_bic(imp, n_trees))

class testGetNewTrees(unittest.TestCase):

    def test1_t2_to_t3(self):
        x = {k:v for k, v in test1.items() if k.startswith("2")}
        n_trees=3
        self.assertEqual(get_new_trees(x, n_trees), [['B', 'AA', 'AB'], ['A', 'BA', 'BB']])

    def test1_t3_to_t4(self):
        """
        As '3_mast_B_AA_AB': None, there should only be 3 new trees instead of 6 
        """
        n_trees=4
        self.assertEqual(
            get_new_trees(test1, n_trees), 
            [['BA', 'BB', 'AA', 'AB'], ['A', 'BB', 'BAA', 'BAB'], ['A', 'BA', 'BBA', 'BBB']]
        )

    def test_truepos_t3_to_t4(self):
        n_trees=4
        self.assertEqual(
            get_new_trees(imp, n_trees), 
            [['AA', 'AB', 'BA', 'BB'], ['B', 'AB', 'AAA', 'AAB'], ['B', 'AA', 'ABA', 'ABB'], 
            ['BA', 'BB', 'AA', 'AB'], ['A', 'BB', 'BAA', 'BAB'], ['A', 'BA', 'BBA', 'BBB']]
        )

class testNewPartitionsPerTree(unittest.TestCase):

    def test_t2_to_t3(self):
        self.assertEqual(new_partitions_per_tree('A'), ['AA', 'AB'])
        self.assertEqual(new_partitions_per_tree('B'), ['BA', 'BB'])

    def test_t3_to_t4(self):
        self.assertEqual(new_partitions_per_tree('AB'), ['ABA', 'ABB'])
        self.assertEqual(new_partitions_per_tree('BA'), ['BAA', 'BAB'])

class testPreviousRun(unittest.TestCase):

    def test_splits(self):
        self.assertEqual(get_previous_run(3, ['A', 'B']), '2_mast_A_B')
        self.assertEqual(get_previous_run(4, ['A', 'BA', 'BB']), '3_mast_A_BA_BB')

class testGetPreviousAln(unittest.TestCase):

    def test_test1(self):
        self.assertEqual(get_previous_aln('2_mast_A_B', 'A', test1), 'test1_class_1-out.fas')

    def test_test1(self):
        self.assertEqual(get_previous_aln('2_mast_A_B', 'A', test1), 'test1_class_1-out.fas')

if __name__ == '__main__':
    unittest.main()
