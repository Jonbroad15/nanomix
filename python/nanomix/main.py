import argparse
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

from functions import *
BISULFITE_ATLAS = os.path.join(script_dir, '..', '..', 'atlases', '39Bisulfite.tsv')

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""
Nanomix: a tool for cell-type deconvolution of mixed methylation calls from Oxford Nanopore sequencing data.
------------------------------------------------------------------------------------------------------------
Version:   v0.1.0
About:     Developed in the simpsonlab.github.io by Jonathan Broadbent.
Code:      https://github.com/simpsonlab/nanomix

This package incorperates 4 different deconvolution models (nnls, llse, llsp, mmse).
Deconvolution results are dependent on the reference atlas used.
We recommend using the atlas from the paper: https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1.full
Which is labelled as '39Bisulfite.tsv' in the 'atlas' directory.

Nanomix sub-commands include:

Deconvolute:    Find the proportion of cell types present in a methylome
Evaluate:       Evaluate the accuracy of a deconvolution model on a given methylome with labelled cell types
Simulate:       Generate a methylome from a given cell-type mixture
Assign:         Assign reads in the methylome to cell-types

All outputs are directed to stdout. Select a sub-command with the [-h] flag for info on hparams
                                     """)
    subparsers = parser.add_subparsers(dest='command', required=True,
                                       help="nanomix functions")
    parser_deconvolute = subparsers.add_parser('deconvolute',formatter_class=argparse.RawDescriptionHelpFormatter,
description="""
Deconvolute a methylome. Return the mixture proportions (sigma) in a tsv format to stdout.

Methylome: A tsv file containing the methylome to be deconvoluted. The file must contain the following columns:
    chr:                The chromosome the read is mapped to
    start:              The start position of the read
    end:                The end position of the read
    total_calls:        The total number of reads mapped to the region
    modified_calls:     The number of modified reads mapped to the region
This file can be created from a BAM file using our associated program, mbtools (https://github.com/jts/mbtools)

Deconvolution Models:
    llse (default):     log-likelihood with sequencing errors. Maximize the likelihood of the methylome, atlas and sigma
                        by assuming modification calls follow a binomial distribution. Good for sequencing data with high error
                        rate and non-uniform coverage distribution. (Oxford Nanopore)
    nnls:               non-negative least squares. Minimize the squared error between the methylome and what we expect for
                        the methylome (given sigma and the atlas). Recommended for fast deconvolution of methylomes with high
                        coverage. (Methylation Arrays)
    mmse:               mixture model with sequencing errors. Also follows a binomial distribution, but softly assigns fragments
                        to cell-types. Optimization uses expectation maximization (slower than above). Recommended for high resolution
                        deconvolution (many cell types) and an atlas with large regions of grouped CpGs.
    llsp:               log-likelihood with sequencing perfect. Same as llse, without error modelling. Useful for differentiating the
                        effect of sequencing errors on deconvolution loss and accuracy.
""")
    parser_deconvolute.add_argument('methylome', help='Path to methylome tsv file with columns: {chr, start, end, total_calls, modified_calls}')
    parser_deconvolute.add_argument('-a', '--atlas', type=str, default=BISULFITE_ATLAS, help='Path to reference atlas')
    parser_deconvolute.add_argument('-p01', default=0.05, type=float, help='Sequencing miscall rate')
    parser_deconvolute.add_argument('-p11', default=0.95, type=float, help='Sequencing correct call rate')
    parser_deconvolute.add_argument('-m', '--model', default='llse', type=str, help='Deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_deconvolute.add_argument('-i', '--sigma_init', default='null', type=str, help='Initalize sigma with one of the other models (mmse only)')
    parser_deconvolute.add_argument('-n', '--max_iter', default=100, type=int, help='Maximum number of iterations for the EM optimization (mmse only)')
    parser_deconvolute.add_argument('-p', '--min_proportion', default=0.01, type=float, help='Minimum proportion of a cell type to be considered (mmse only)')
    parser_deconvolute.add_argument('-t', '--stop_thresh', default=1e-3, type=float, help='Stop EM iterations when percent change in log-likelihood falls below this value (mmse only)')
    parser_deconvolute.set_defaults(func=deconvolute)

    parser_evaluate = subparsers.add_parser('evaluate', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Evaluate the loss and accuracy of cell-type deconvolution on a given methylome with known cell types.

Metrics outputed to stdout:
    Deconvolution loss:         Euclidean distance between the computed and true mixture proportion.
    Accuracy at confidence <t>: Proportion of correctly assigned reads with confience >= t.


Methylome: A tsv file containing the methylome to be deconvoluted. The file must contain the following columns:
    chr:                The chromosome the read is mapped to
    start:              The start position of the read
    end:                The end position of the read
    total_calls:        The total number of reads mapped to the region
    modified_calls:     The number of modified reads mapped to the region
    cell_type:          String identifier of cell_type that matches that in the atlas header
The simulate function can be used to create this file.
""")
    parser_evaluate.add_argument('methylome', help='Path to methylome tsv file with columns: {chr, start, end, total_calls, modified_calls, cell_type}')
    parser_evaluate.add_argument('-a', '--atlas', type=str, default=BISULFITE_ATLAS, help='Path to reference atlas')
    parser_evaluate.add_argument('-p01', default=0.05, type=float, help='Sequencing miscall rate')
    parser_evaluate.add_argument('-p11', default=0.95, type=float, help='Sequencing correct call rate')
    parser_evaluate.add_argument('-m', '--model', default='llse', type=str, help='Deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_evaluate.set_defaults(func=evaluate)

    parser_simulate = subparsers.add_parser('simulate', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Simulate a methylome from a given cell-type mixture (sigma).

Sigma: A tsv file containing the cell-type mixture proportions. The file must contain the following columns:
    cell_type:          String identifier of cell_type that matches that in the atlas header
    proportion:         float between 0 and 1
""")
    parser_simulate.add_argument('sigma', type=str, help='Path to sigma tsv file')
    parser_simulate.add_argument('-a', '--atlas', type=str, default=BISULFITE_ATLAS, help='Path to reference atlas')
    parser_simulate.add_argument('-p01', default=0.05, type=float, help='Sequencing miscall rate')
    parser_simulate.add_argument('-p11', default=0.95, type=float, help='Sequencing correct call rate')
    parser_simulate.add_argument('-c', '--coverage', default=1, type=float, help='Sequencing coverage')
    parser_simulate.add_argument('-r', '--region_size', default=5, type=int, help='Number of CpGs in each region')
    parser_simulate.set_defaults(func=simulate)

    parser_assign = subparsers.add_parser('assign', formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Assign cell types to a methylome using a given atlas and pre-deconvoluted sigma.
Output (to stdout) is of the same format as the methylome with an appended column indicating the cell_type.

Methylome:  A tsv file containing the methylome to be deconvoluted. The file must contain the following columns:
    chr:                The chromosome the read is mapped to
    start:              The start position of the read
    end:                The end position of the read
    total_calls:        The total number of reads mapped to the region
    modified_calls:     The number of modified reads mapped to the region
    cell_type:          String identifier of cell_type that matches that in the atlas header
The simulate function can be used to create this file.
""")
    parser_assign.add_argument('methylome', help='Path to methylome tsv file with columns: {chr, start, end, total_calls, modified_calls}')
    parser_assign.add_argument('-a', '--atlas', type=str, default=BISULFITE_ATLAS, help='Path to reference atlas')
    parser_assign.add_argument('-p01', default=0.05, type=float, help='Sequencing miscall rate')
    parser_assign.add_argument('-p11', default=0.95, type=float, help='Sequencing correct call rate')
    parser_assign.add_argument('-s', '--sigma', required=True, type=str, help='Path to sigma tsv file')
    parser_assign.set_defaults(func=assign_fragments)

    args = parser.parse_args()
    arg_dict = {k : v for k, v in vars(args).items() if k not in ['func','command']}
    args.func(**arg_dict)

if __name__ == "__main__":
    main()
