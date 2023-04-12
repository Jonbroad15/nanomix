# Test functions in main.py
import unittest
import os
import sys


script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

from main import *

script_dir = os.path.dirname(os.path.realpath(__file__))
methylome = os.path.join(script_dir, 'test_data', 'test_methylome.tsv')
atlas = os.path.join(script_dir, 'test_data', 'test_atlas.tsv')
sigma_path = os.path.join(script_dir, 'test_data', 'test_sigma.tsv')

class TestMain(unittest.TestCase):

    @classmethod
    def setupclass(cls):
        # write sigma
        sigma = [0.3, 0.7]
        coverage = 100
        # write an atlas
        with open(atlas, 'w') as f:
            f.write('chr\tstart\tend')
            for i in range(len(sigma)):
                f.write(f'\tcell{i}')
            f.write('\n')
            for i in range(len(sigma)):
                f.write(f'chr1\t{i*100}\t{i*100+99}')
                for j in range(len(sigma)):
                    if i == j:
                        f.write('\t0.99')
                    else:
                        f.write('\t0.01')
                f.write('\n')
        # simulate a methylome
        with open(methylome, 'w') as f:
            f.write('chr\tstart\tend\ttotal_calls\tmodified_calls\n')
            for i in range(len(sigma)):
                # simulate n reads for every region
                for _ in range(coverage):
                    # choose a cell type wighted by sigma
                    cell_type = np.random.choice(len(sigma), p=sigma)
                    if cell_type == i:
                        f.write(f'chr1\t{i*100}\t{i*100+99}\t1\t1\n')
                    else:
                        f.write(f'chr1\t{i*100}\t{i*100+99}\t1\t0\n')

    # def test_llse(self):
        # # test deconvolute function
        # cell_type_proportions = deconvolute(methylome, atlas, 'llse')

        # # check that the cell_type proportion values are close
        # self.assertAlmostEqual(cell_type_proportions['cell1'], 0.7, places=1)

    # def test_nnls(self):
        # cell_type_proportions = deconvolute(methylome, atlas, 'nnls')

        # # check that the cell_type proportion values are close
        # self.assertAlmostEqual(cell_type_proportions['cell1'], 0.7, places=1)

    # def test_llsp(self):
        # cell_type_proportions = deconvolute(methylome, atlas, 'llsp')

        # # check that the cell_type proportion values are close
        # self.assertAlmostEqual(cell_type_proportions['cell1'], 0.7, places=1)

    def test_mmse(self):
        cell_type_proportions = deconvolute(methylome, atlas, 'mmse')

        # check that the cell_type proportion values are close
        self.assertAlmostEqual(cell_type_proportions['cell1'], 0.7, places=1)



if __name__ == '__main__':
    unittest.main()
