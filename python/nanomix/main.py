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
    :return: None
    """
    #TODO: Add evaluation code

def fit_model(model, atlas, sample, p01, p11):
    """
    Wrapper function to select model for deconvolution

    :param model: model to fit
    :param atlas: reference atlas
    :param sample: sample to deconvolute
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: None
    """
    if model == 'nnls':
        sigma = fit_nnls(atlas, s)
    elif model == 'llse':
        sigma = fit_llse(atlas, s, p01, p11)
    elif model == 'llsp':
        sigma = fit_llsp(atlas, s, p01, p11)
    elif model == 'null':
        sigma = fit_uniform(atlas, s)
    else:
        raise ValueError(f"no such model: {model}. Choose from [nnls, llse, llsp, mmse]")
    return sigma

def deconvolute(methylomes, atlas, model, p01, p11, sigma_init, max_iter, min_proportion, stop_thresh):
    """
    Deconvolute a methylome using a given model to get the proportion of each cell type present

    :param methylome: Path to tsv file of methylome
    :param atlas: Path to tsv file of atlas
    :param model: Deconvolution model options: [nnls, llse, llsp, mmse]
    :param p01: Sequencing miscall rate
    :param p11: Sequencing correct call rate
    :param sigma_init: Initalize sigma with one of the other models (mmse only)
    :param max_iter: Maximum number of iterations for the model (mmse only)
    :param min_proportion: Minimum proportion of a cell type to be considered (mmse only)
    :param stop_thresh: Threshold for stopping iterations (mmse only)
    :return: none
    """

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
        Exception("Empty methylome file")
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
    if model == 'mmse':
        #TODO: change output of fit_mmse to be consistent with other models (dictionary)
        sigma = fit_model(init, atlas, s, p01, p11)
        cell_type_proportions = fit_mmse(atlas, methylome, sigma, p01, p11, stop_thresh=1e-3, max_iter=100, min_proportion=0.01)
    else:
        sigma = fit_model(model, atlas, s, p01, p11)
        cell_type_proportions = {cell_type: proportion for cell_type, proportion in zip(atlas.get_cell_types(), sigma)}

    # return deconvolution results
    print("cell_type\tproportion")
    for cell_type, proportion in cell_type_proportions.items():
        print(f"{cell_type}\t{proportion}")

    # TODO: return cell type assignments



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--atlas', type=str, default='atlases/39Bisulfite.csv')
    parser.add_argument('-p01', default=0.05, type=float, help='Sequencing miscall rate')
    parser.add_argument('-p11', default=0.95, type=float, help='Sequencing correct call rate')
    parser.add_argument('methylome', help='Path to methylome tsv file')
    subparsers = parser.add_subparsers(dest='command', required=True)

    parser_deconvolute = subparsers.add_parser('deconvolute')
    parser_deconvolute.add_argument('-m', '--model', default='llse', type=str, help='Deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_deconvolute.add_argument('-i', '--sigma_init', default='null', type=str, help='Initalize sigma with one of the other models (mmse only)')
    parser_deconvolute.add_argument('-n', '--max_iter', default=100, type=int, help='Maximum number of iterations for the model (mmse only)')
    parser_deconvolute.add_argument('-p', '--min_proportion', default=0.01, type=float, help='Minimum proportion of a cell type to be considered (mmse only)')
    parser_deconvolute.add_argument('-t', '--stop_thresh', default=1e-3, type=float, help='Stop EM iterations when change in log-likelihood falls below this value (mmse only)')
    parser_deconvolute.set_defaults(func=deconvolute)

    parser_evaluate = subparsers.add_parser('evaluate')
    parser_evaluate.add_argument('-m', '--model', default='nnls', type=str, help='Deconvolution model options: [nnls, llse, llsp, mmse]')
    parser_evaluate.set_defaults(func=evaluate)

    parser_simulate = subparsers.add_parser('simulate')
    parser_simulate.add_argument('-c' '--coverage', default=1, type=float, help='Sequencing coverage')
    parser_simulate.add_argument('-r', '--region_size', default=5, help='Number of CpGs in each region')
    parser.add_argument('-s', '--sigma', required=True, type=str, help='Path to sigma tsv file')
    parser_simulate.set_defaults(func=simulate)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
