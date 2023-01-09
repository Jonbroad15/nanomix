#! /usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import pyranges as pr

from nanomix import nanomix
from models import fit_llse, fit_llsp, fit_nnls, fit_mmse, fit_uniform
from atlas import ReferenceAtlas, Sample

def simulate(methylome, atlas, sigma, coverage, region_size, p01, p11):
    """
    Simulate a methylome from a reference atlas at the given cell type proportion
    Wrapper function to call Rust code from Python

    :param methylome: path to save simulated methylome
    :param atlas: path to reference atlas
    :param sigma: path to tsv file containing cell type proportions
    :param coverage: coverage of simulated reads
    :param region_size: number of CpGs in each region
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: None
    """
    nanomix.generate_methylome(methylome, atlas, sigma, coverage, region_size, p01, p11)

def evaluate(methylome, atlas, model, p01, p11):
    """
    Evaluate the performance of a model on a simulated methylome
    Wrapper function to call Rust code from Python

    :param methylome: path to simulated methylome
    :param atlas: path to reference atlas
    :param model: model to evaluate
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    """
    #TODO: Add evaluation code

def deconvolute(methylomes, atlas, model, p01, p11):
    """
    Deconvolute a methylome using a given model to get the proportion of each cell type present

    :param methylome: path to tsv file of methylome
    :param atlas: path to tsv file of atlas
    :param model: deconvolution model options: [nnls, llse, llsp, mmse]
    :param p01: sequencing miscall rate
    :param p11: sequencing correct call rate
    :param random_inits: use random initializations for llse
    :return: deconvolution results
    """

    if model == 'mmse':
        #TODO: add initializations for mmse
        #TODO: add extra params for MMSE
        cell_type_proportions = fit_mmse(atlas, methylome, p01, p11, stop_thresh=1e-3, max_iter=100, min_proportion=0.01)
    else:
        # load atlas
        columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                                'start':'Start',
                                'end':'End'}
        df_atlas = pd.read_csv(atlas, sep='\t').rename(columns=columns)
        df_atlas.drop_duplicates(inplace=True)
        if 'label' in df_atlas.columns: df_atlas.drop('label', axis=1, inplace=True)
        df_atlas.dropna(inplace=True)
        gr_atlas = pr.PyRanges(df_atlas).sort()

        # Read methylomes data from mbtools
        try:
            df = pd.read_csv(methylome, sep='\t').rename(columns=columns)
        except pd.errors.EmptyDataError:
            Exception("Empty Methylome file")
        df.dropna(inplace=True)
        gr_sample = pr.PyRanges(df).sort()

        # Join atlas and sample
        gr = gr_atlas.join(gr_sample)
        atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])
        t = np.array(gr.total_calls, dtype=np.float32)
        m = np.array(gr.modified_calls, dtype=np.float32)

        xhat = m/t
        name = get_sample_name(methylome)
        sample_names.append(name)
        s = Sample(name, xhat, m, t)

        # Run
        if model == 'nnls':
            sigma = fit_nnls(atlas, s)
        elif model == 'llse':
            sigma = fit_llse(atlas, s, p01, p11, random_inits)
        elif model == 'llsp':
            sigma = fit_llsp(atlas, s, p01, p11, random_inits)
        elif model == 'null':
            sigma = fit_uniform(atlas, s)
        else:
            Exception(f"no such model {model}")

        cell_type_proportions = {cell_type: proportion for cell_type, proportion in zip(atlas.get_cell_types(), sigma)}
    # TODO: Print output of cell type proportions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--atlas', type=str, default='atlases/39Bisulfite.csv')
    parser.add_argument('-p01', default=0.05, type=float, help='sequencing miscall rate')
    parser.add_argument('-p11', default=0.95, type=float, help='sequencing correct call rate')
    parser.add_argument('methylome', help='path to methylome tsv file')
    subparsers = parser.add_subparsers(dest='command', required=True)

    parser_deconvolute = subparsers.add_parser('deconvolute')
    parser_deconvolute.add_argument('-m', '--model', default='nnls', type=str, help='deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_deconvolute.set_defaults(func=deconvolute)

    parser_evaluate = subparsers.add_parser('evaluate')
    parser_evaluate.add_argument('-m', '--model', default='nnls', type=str, help='deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_evaluate.set_defaults(func=evaluate)

    parser_simulate = subparsers.add_parser('simulate')
    parser_simulate.add_argument('-c' '--coverage', default=1, type=float, help='sequencing coverage')
    parser_simulate.add_argument('-r', '--region_size', default=5, help='number of CpGs in each region')
    parser.add_argument('-s', '--sigma', required=True, type=str, help='path to sigma tsv file')
    parser_simulate.set_defaults(func=simulate)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
