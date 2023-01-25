# Test functions in main.py
import unittest
import os
import sys


script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

from main import *
from tools import *

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
        cell_types = get_cell_types(atlas)
        with open(methylome, 'w') as f:
            f.write('chr\tstart\tend\ttotal_calls\tmodified_calls\tcell_type\n')
            for i in range(len(sigma)):
                # simulate n reads for every region
                for _ in range(coverage):
                    # choose a cell type wighted by sigma
                    cell_type = np.random.choice(len(sigma), p=sigma)
                    if cell_type == i:
                        f.write(f'chr1\t{i*100}\t{i*100+99}\t1\t1\t{cell_types[cell_type]}\n')
                    else:
                        f.write(f'chr1\t{i*100}\t{i*100+99}\t1\t0\t{cell_types[cell_type]}\n')

    def test_evaluate(self):
        # Test evaluate function
        evaluate(methylome, atlas, 'mmse')

if __name__ == '__main__':
    unittest.main()
