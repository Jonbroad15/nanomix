# Test the joining operation on the atlas and methylome data
import unittest
import os
import sys


script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

from tools import *

script_dir = os.path.dirname(os.path.realpath(__file__))
methylome = os.path.join(script_dir, 'test_data', 'test_methylome.tsv')
atlas = os.path.join(script_dir, 'test_data', 'test_atlas.tsv')
sigma_path = os.path.join(script_dir, 'test_data', 'test_sigma.tsv')

class TestMain(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        sigma = [0.3, 0.7]
        coverage = 3
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

    def test_load_atlas(self):
        reference_atlas, sample = load_atlas(atlas, methylome)
        self.assertEqual(len(sample.m), 6)
        self.assertEqual(len(reference_atlas.A), 6)
        print(reference_atlas.A)




if __name__ == '__main__':
    unittest.main()
