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
    def setUpClass(cls):
        # Write sigma
        sigma = [0.3, 0.7]
        coverage = 100
        with open(sigma_path, 'w') as f:
            f.write('cell_type\tproportion\n')
            for i, s in enumerate(sigma):
                f.write(f'cell{i}\t{s}\n')
        # Write an atlas
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
        # Simulate a methylome
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

    def test_simulate(self):
        # Test simulate function
        simulate(atlas, sigma_path, 200, 5, 0.0, 1.0)

        # check that the cell_type proportions are close to sigma
        # Tissues are enumerated as 0 and 1. so the average of them should be the proportion of tissue2
        # with open(methylome, 'r') as f:
            # lines = f.readlines()
            # cell_types = [int(line.split('\t')[-1].strip()) for line in lines[1:]]
        # self.assertAlmostEqual(sum(cell_types)/len(cell_types), 0.7, places=1)

if __name__ == '__main__':
    unittest.main()
