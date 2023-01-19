#! /usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import pyranges as pr
import math
from scipy.stats import binom, dirichlet
from scipy.optimize import minimize, nnls, Bounds
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

from _nanomix import *
from atlas import ReferenceAtlas, Sample
from tools import *


def fit_model(methylome, atlas_path, model, p01, p11):
    """
    Wrapper function to select model for deconvolution
    mmse model is initalized seperately with Rust bindings
    methylome and atlas are joined on covered regions before deconvolution

    :param methylome: path to simulated methylome
    :param atlas_path: path to reference atlas
    :param model: model to fit
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: None
    """
    # load atlas
    columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                            'start':'Start',
                            'end':'End'}
    df_atlas = pd.read_csv(atlas_path, sep='\t').rename(columns=columns)
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
    # Check for empty upon join
    if len(gr) == 0:
        Exception("Empty join between atlas and sample. The sample does not overlap with any regions in the atlas.")
    atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])
    t = np.array(gr.total_calls, dtype=np.float32)
    m = np.array(gr.modified_calls, dtype=np.float32)

    xhat = m/t
    s = Sample('methylome', xhat, m, t)
    if model == 'nnls':
        sigma = fit_nnls(atlas, s)
    elif model == 'llse':
        sigma = fit_llse(atlas, s, p01, p11)
    elif model == 'llsp':
        sigma = fit_llsp(atlas, s)
    elif model == 'null':
        sigma = fit_uniform(atlas, s)
    else:
        raise ValueError(f"no such model: {model}. Choose from [nnls, llse, llsp, mmse]")
    return sigma

def log_likelihood_sequencing_perfect(atlas, sigma, sample):
    """
    Compute the log-likelihood of the llsp model

    :param atlas: reference atlas
    :param sigma: cell-type proportions
    :param sample: sample methylome
    :return: log-likelihood
    :rtype: float
    """
    sigma_t = sigma.reshape( (atlas.K, 1) )
    x = np.clip(np.ravel(atlas.get_x(sigma_t)), 0, 1.0)
    b =  binom.logpmf(sample.m, sample.t, x)
    return np.sum(b)

def fit_llsp(atlas, sample,
             n_trials=10):
    """
    fit the log-likelihood sequencing perfect model

    :param atlas: reference atlas
    :param sample: sample methylome
    :param n_trials: number of trials to run with random initalizations of sigma
    :return: cell-type proportions
    :rtype: np.array
    """
    f = lambda x: -1 * log_likelihood_sequencing_perfect(atlas, x, sample)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / atlas.K ] * atlas.K)
    n_trials = 10
    best_ll = np.inf
    best_sol = None
    initializations = dirichlet.rvs(alpha, size=n_trials).tolist()

    for (i, init) in enumerate(initializations):
        res = minimize(f, init, method='SLSQP', options={'maxiter': 100, 'disp':False}, bounds=bnds, constraints=cons)
        ll = res.get("fun")
        if ll < best_ll:
            best_ll = ll
            best_sol = res
    return best_sol.x/np.sum(best_sol.x)

def log_likelihood_sequencing_with_errors(atlas, sigma, sample, p01,p11):
    """
    Compute the log-likelihood of the llse model

    :param atlas: reference atlas
    :param sigma: cell-type proportions
    :param sample: sample methylome
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: log-likelihood
    :rtype: float
    """
    sigma_t = sigma.reshape( (atlas.K, 1) )

    # the solver we use can try values that are outside
    # the constraints we impose, we need to clip here to prevent
    # things from blowing up
    x = np.clip(np.ravel(atlas.get_x(sigma_t)), 0, 1.0)
    if p11:
        p = x * p11 + (1-x) * p01
    else:
        p = x * (1 - p01) + (1 - x) * p01
    b =  binom.logpmf(sample.m, sample.t, p)
    binomial_coef = sum([math.log(math.comb(int(t), int(m))) for m,t in zip(sample.m, sample.t)])

    return np.sum(b) - binomial_coef

def fit_uniform(atlas, sample):
    """
    fit the uniform model (null model)

    :param atlas: reference atlas
    :param sample: sample methylome
    :return: cell-type proportions
    :rtype: np.array
    """
    return np.array([1.0 / atlas.K ] * atlas.K)

def fit_llse(atlas, sample, p01, p11,
             n_trials=10):
    """
    fit the log-likelihood sequencing with errors model

    :param atlas: reference atlas
    :param sample: sample methylome
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :param n_trials: number of trials to run with random initalizations of sigma
    :return: cell-type proportions
    :rtype: np.array
    """
    f = lambda x: -1 * log_likelihood_sequencing_with_errors(atlas, x, sample, p01, p11)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / atlas.K ] * atlas.K)
    best_ll = np.inf
    best_sol = None
    initializations = dirichlet.rvs(alpha, size=n_trials).tolist()

    for (i, init) in enumerate(initializations):
        res = minimize(f, init, method='SLSQP', options={'maxiter': 100, 'disp':False}, bounds=bnds, constraints=cons)
        ll = res.get("fun")
        if ll < best_ll:
            best_ll = ll
            best_sol = res
    return best_sol.x/np.sum(best_sol.x)

def fit_nnls(atlas, sample):
    """
    fit the non-negative least squares model

    :param atlas: reference atlas
    :param sample: sample methylome
    :return: cell-type proportions
    :rtype: np.array
    """

    # add sum=1 constraint
    t = np.array([1.0] * atlas.K).reshape( (1, atlas.K) )
    A = np.append(atlas.A, t, axis=0)
    b = np.append(sample.x_hat, [1.0], axis=0)
    res = nnls(A, b)
    return res[0]/np.sum(res[0])

def fit_mmse(atlas, sample, sigma, p01, p11, stop_thresh, max_iter, min_proportion,
             true_sigma=None, true_assignments=None):
    """
    Fit mixture model with sequencing errors (MMSE) onto sample with reference atlas.
    Wrapper function to access Rust implementation.

    :param atlas: path to reference atlas
    :param sample: path to sample methylome to fit
    :param p01: sequencing miscall rate
    :param p11: sequencing correct call rate
    :param stop_thresh: threshold for stopping iterations
    :param max_iter: maximum number of iterations
    :param min_proportion: minimum proportion of a cell type in the mixture
    :return: fitted mixture model
    """
    # initialize MMSE model
    mmse = MMSE(sample, atlas, sigma, p01, p11)
    if true_sigma is not None and true_assignments is not None:
        mmse.evaluate(stop_thresh, max_iter, min_proportion, true_sigma, true_assignments)
    else:
        mmse.optimize(stop_thresh, max_iter, min_proportion)
    return mmse.cell_type_proportions()

