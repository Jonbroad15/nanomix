#! /usr/bin/env python

import argparse
import numpy as np
import sys
import csv
import os
import re
import pandas as pd
import pyranges as pr
import math
from nanomix import nanomix

from scipy.stats import binom, dirichlet
from scipy.optimize import minimize, nnls, Bounds

script_dir = os.path.dirname(__file__)
ATLAS = os.path.join(script_dir, '..', 'atlases', 'meth_atlas.csv')

class ReferenceAtlas:
    def __init__(self, gr):
        self.cpg_ids = [(chrom, start, end) for chrom, start, end in\
                        zip(gr.Chromosome, gr.Start, gr.End)]
        cell_types = list(gr.columns[3:])
        self.K = len(cell_types)
        self.v = {k:list(gr[k]) for k in cell_types}
        self.A = np.array(gr.loc[:, list(cell_types)])

    def get_x(self, sigma):
        x = np.matmul(self.A, sigma)
        return x

    def get_num_cpgs(self):
        return len(self.cpg_ids)

    def get_cell_types(self):
        return list(self.v.keys())

    def get_num_cell_types(self):
        return len(self.v.keys())

class Sample:
    def __init__(self, name, x_hat, m, t):
        self.name = name
        self.x_hat = x_hat
        self.m = m
        self.t = t

def eq_constraint(x):
    return 1 - np.sum(x)

#
# Model wrappers
#
def log_likelihood_sequencing_perfect(atlas, sigma, sample, p01,p11):
    sigma_t = sigma.reshape( (atlas.K, 1) )
    x = np.clip(np.ravel(atlas.get_x(sigma_t)), 0, 1.0)
    b =  binom.logpmf(sample.m, sample.t, x)
    return np.sum(b)

def eq_constraint(x):
    return 1 - np.sum(x)
def fit_llsp(atlas, sample, p01, p11, random_inits):
    f = lambda x: -1 * log_likelihood_sequencing_perfect(atlas, x, sample, p01, p11)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / atlas.K ] * atlas.K)
    if random_inits:
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
    else:
        res = minimize(f, alpha, method='SLSQP', options={'maxiter': 100, 'disp':False}, bounds=bnds, constraints=cons)
        return res.x / np.sum(res.x)

# Binomial model with sequencing errors, when p01 = 0
# this is the same as the perfect data model
def log_likelihood_sequencing_with_errors(atlas, sigma, sample, p01,p11):
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
    return np.array([1.0 / atlas.K ] * atlas.K)

def fit_llse(atlas, sample, p01, p11, random_inits):
    f = lambda x: -1 * log_likelihood_sequencing_with_errors(atlas, x, sample, p01, p11)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / atlas.K ] * atlas.K)
    if random_inits:
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
    else:
        res = minimize(f, alpha, method='SLSQP', options={'maxiter': 100, 'disp':False}, bounds=bnds, constraints=cons)
        return res.x / np.sum(res.x)

def fit_nnls(atlas, sample):

    # add sum=1 constraint
    t = np.array([1.0] * atlas.K).reshape( (1, atlas.K) )
    A = np.append(atlas.A, t, axis=0)
    b = np.append(sample.x_hat, [1.0], axis=0)
    res = nnls(A, b)
    return res[0]/np.sum(res[0])

def fit_nnls_constrained(atlas, sample):
    sigma_0 = np.array([ [ 1.0 / atlas.K ] * atlas.K ])
    f = lambda x: np.linalg.norm(atlas.A.dot(x) - sample.x_hat)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    res = minimize(f, sigma_0, method='SLSQP', options={'maxiter': 10, 'disp':False}, bounds=bnds, constraints=cons)
    return res.x

def fit_mmse(atlas, sample, p01, p11, stop_thresh, max_iter, min_proportion):
    """
    Fit mixture model with sequencing errors (MMSE) onto sample with reference atlas

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
    mmse = nanomix.MMSE(sample, atlas, p01, p11)
    mmse.optimize(stop_thresh, max_iter, min_proportion)

    return mmse.sigma

