#! /usr/bin/env python

import argparse
import numpy as np
import sys
import csv
import os
import pandas as pd
import pyranges as pr

from scipy.stats import binom, dirichlet
from scipy.optimize import minimize, nnls, Bounds

script_dir = os.path.dirname(__file__)
ATLAS = os.path.join(script_dir, '..', 'atlases', 'meth_atlas.csv')

class ReferenceAtlas:
    def __init__(self, gr):
        self.cpg_ids = [(chrom, start, end) for chrom, start, end in\
                        zip(gr.Chromosome, gr.Start, gr.End)]
        cell_types = set(gr.columns) - {'Chromosome', 'Start', 'End', 'type'}
        self.K = len(cell_types)
        self.v = {k:list(gr[k]) for k in cell_types}
        self.A = np.array(gr.loc[:, cell_types])

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

# Binomial model with sequencing errors, when epsilon = 0
# this is the same as the perfect data model
def log_likelihood_sequencing_with_errors(atlas, sigma, sample, epsilon):
    sigma_t = sigma.reshape( (atlas.K, 1) )

    # the solver we use can try values that are outside
    # the constraints we impose, we need to clip here to prevent
    # things from blowing up
    x = np.clip(np.ravel(atlas.get_x(sigma_t)), 0, 1.0)
    p = x * (1 - epsilon) + (1 - x) * epsilon
    b =  binom.logpmf(sample.m, sample.t, p)

    #print("SigmaT", sigma_t)
    #print("SigmaSum", np.sum(sigma_t))
    #print("m", sample.m)
    #print("t", sample.t)
    #print("x", x)
    #print("B", b)
    #print("Sum", np.sum(b))
    return np.sum(b)

def eq_constraint(x):
    return 1 - np.sum(x)

#
# Model wrappers
#
def fit_llse(atlas, sample, epsilon):
    n_trials = 10

    f = lambda x: -1 * log_likelihood_sequencing_with_errors(atlas, x, sample, epsilon)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / atlas.K ] * atlas.K)
    best_ll = np.inf
    best_sol = None

    #initializations = dirichlet.rvs(alpha, size=n_trials).tolist()
    initializations = [ alpha ] # uniform

    for (i, init) in enumerate(initializations):
        sigma_0 = dirichlet.rvs(alpha, size=1)
        res = minimize(f, init, method='SLSQP', options={'maxiter': 10, 'disp':False}, bounds=bnds, constraints=cons)
        ll = res.get("fun")
        if ll < best_ll:
            best_ll = ll
            best_sol = res
    return best_sol.x

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

def get_sample_name(s):
    s = s.split('/')[-1]
    s = s.split('.')[0]
    return s

def deconvolve(methylomes, atlas, model, epsilon):
    Y = []
    sample_names = []
    columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                            'start':'Start',
                            'end':'End'}
    df_atlas = pd.read_csv(atlas, sep='\t').rename(columns=columns)
    df_atlas.drop_duplicates(inplace=True)
    df_atlas.dropna(inplace=True)
    gr_atlas = pr.PyRanges(df_atlas).sort()
    for methylome in methylomes:
        # read methylomes data from mbtools
        try:
            df = pd.read_csv(methylome, sep='\t').rename(columns=columns)
        except pd.errors.EmptyDataError:
            continue
        df.drop_duplicates(inplace=True)
        df.dropna(inplace=True)
        gr_sample = pr.PyRanges(df).sort()

        # Init atlas and sample
        gr = gr_atlas.join(gr_sample)
        atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])
        xhat = np.array(gr.modification_frequency)
        t = np.array(gr.num_called_reads)
        m = np.rint((t * xhat))
        name = get_sample_name(methylome)
        sample_names.append(name)
        s = Sample(name, xhat, m, t)

        # Run
        if model == 'nnls':
            Y.append(fit_nnls(atlas, s))
        else:
            Y.append(fit_llse(atlas, s, epsilon))

    return Y, sample_names, atlas

def deconvolve_uxm(methylomes, atlas, model, epsilon):
    Y = []
    sample_names = []
    columns={'chromosome':'Chromosome',
                            'start':'Start',
                            'end':'End'}
    df_atlas = pd.read_csv(atlas, sep='\t').rename(columns=columns)
    df_atlas.drop_duplicates(inplace=True)
    gr_atlas = pr.PyRanges(df_atlas).sort()
    for methylome in methylomes:
        # read methylomes data from mbtools
        try:
            df = pd.read_csv(methylome, sep='\t').rename(columns=columns)
        except pd.errors.EmptyDataError:
            continue
        df.drop_duplicates(inplace=True)
        df.dropna(inplace=True)
        gr_sample = pr.PyRanges(df).sort()

        # Init atlas and sample
        gr = gr_atlas.join(gr_sample)
        atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])

        # Combine U/M reads into methylome
        breakpoint()
        u_reads = (np.array(gr.type) == 'U')*np.array(gr.u_reads)
        m_reads = (np.array(gr.type) == 'M')*np.array(gr.m_reads)
        um_reads = u_reads + m_reads
        # sanity check
        if any(um_reads > m_reads): Exception("Error in U/M read masking")

        # Init sample
        t = np.array(gr.num_called_reads)
        name = get_sample_name(methylome)
        sample_names.append(name)
        s = Sample(name, um_reads, um_reads, t)

        # Run
        if model == 'nnls':
            # Scale atlas by coverage in each position
            atlas.A = np.dot(np.diag(t), atlas.A)
            Y.append(fit_nnls(atlas, s))
        else:
            Y.append(fit_llse(atlas, s, epsilon))

    return Y, sample_names, atlas

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', type=str,
            default='/.mounts/labs/simpsonlab/users/jbroadbent/code/cfdna/nanopore_cfdna/atlases/meth_atlas.csv')
    parser.add_argument('--name', type=str, default='sample1')
    parser.add_argument('--model', default='llse', type=str, help='deconvolution model options: [nnml, llse]')
    parser.add_argument('input', nargs='+',
                        help='reference_modifications.tsv file')
    parser.add_argument('--epsilon', default=0.05, type=float)
    parser.add_argument('--uxm', action='store_true', help='Loyfer UXM deconvolution method')
    args = parser.parse_args()

    if args.uxm:
        Y, sample_names, atlas = deconvolve_uxm(args.input, args.atlas, args.model, args.epsilon)
    else:
        Y, sample_names, atlas = deconvolve(args.input, args.atlas, args.model, args.epsilon)

    print("\t".join(['ct'] + sample_names))
    for i, cell_type in enumerate(atlas.get_cell_types()):
        print("\t".join([cell_type] + [str(round(y[i],4)) for y in Y]))

if __name__ == "__main__":
    main()
