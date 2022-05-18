#! /usr/bin/env python

import argparse
import numpy as np
import sys
import csv
import os
import pandas as pd
import pyranges as pr

from scipy.stats import binom
from scipy.optimize import minimize, nnls, Bounds

script_dir = os.path.dirname(__file__)
ATLAS = os.path.join(script_dir, '..', 'atlases', 'meth_atlas.csv')

class ReferenceAtlas:
    def __init__(self, gr):
        self.cpg_ids = [(chrom, start, end) for chrom, start, end in\
                        zip(gr.Chromosome, gr.Start, gr.End)]
        cell_types = list(gr.columns)[3:]
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
    sigma_0 = np.array([ [ 1.0 / atlas.K ] * atlas.K ])
    f = lambda x: -1 * log_likelihood_sequencing_with_errors(atlas, x, sample, epsilon)

    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    res = minimize(f, sigma_0, method='SLSQP', options={'maxiter': 10, 'disp':False}, bounds=bnds, constraints=cons)
    return res.x

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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', type=str,
            default='/.mounts/labs/simpsonlab/users/jbroadbent/code/cfdna/nanopore_cfdna/atlases/meth_atlas.csv')
    parser.add_argument('--name', type=str, default='sample1')
    parser.add_argument('--model', default='llse', type=str, help='deconvolution model options: [nnml, llse]')
    parser.add_argument('input', nargs='+',
                        help='reference_modifications.tsv file')
    parser.add_argument('-o')
    parser.add_argument('--epsilon', default=0.05, type=float)
    parser.add_argument('--fill', action='store_true')
    args = parser.parse_args()

    Y = []
    sample_name = []
    columns={'chromosome':'Chromosome',
                            'start':'Start',
                            'end':'End'}
    df_atlas = pd.read_csv(args.atlas, sep='\t').rename(columns=columns)
    df_atlas.drop_duplicates(inplace=True)
    gr_atlas = pr.PyRanges(df_atlas).sort()
    for input_file in args.input:
        # read input data from mbtools
        try:
            df = pd.read_csv(input_file, sep='\t').rename(columns=columns)
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
        sample_name.append(get_sample_name(input_file))
        s = Sample(args.name, xhat, m, t)

        # Run
        if args.model == 'nnls':
            Y.append(fit_nnls(atlas, s))
        else:
            Y.append(fit_llse(atlas, s, args.epsilon))
    # output
    print("\t".join(['ct'] + sample_name))
    for i, cell_type in enumerate(atlas.get_cell_types()):
        print("\t".join([cell_type] + [str(round(y[i],4)) for y in Y]))

if __name__ == "__main__":
    main()
