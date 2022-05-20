#! /usr/bin/env python

import argparse
import numpy as np
import sys
import csv
import os
import pandas as pd

from scipy.stats import binom, dirichlet
from scipy.optimize import minimize, nnls, Bounds

script_dir = os.path.dirname(__file__)
ATLAS = os.path.join(script_dir, '..', 'atlases', 'meth_atlas.csv')

class ReferenceAtlas:
    def __init__(self, filename, covered_positions):
        self.cpg_ids = list()
        self.v = dict()
        self.K = None
        self.covered_positions = covered_positions

        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            data = list()
            for row in reader:
                chrom = row['chr']
                start = int(row['start'])
                end = int(row['end'])
                cpg_id = (chrom, start, end)
                if cpg_id not in self.covered_positions: continue
                if cpg_id in self.cpg_ids: continue
                self.cpg_ids.append((chrom, start, end))
                cell_types = list(row.keys())[3:]
                self.K = len(cell_types)
                r = list()
                for k in cell_types:
                    if k not in self.v:
                        self.v[k] = list()
                    self.v[k].append(float(row[k]))
                    r.append(float(row[k]))
                data.append(r)

        self.A = np.array(data).reshape((self.get_num_cpgs(), self.get_num_cell_types()))

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
    #s = s.split('.')[1]
    return s
def fill_forward(x):
    prev = 0.0
    for i in range(len(x)):
        if np.isnan(x[i]):
            x[i] = prev
        prev = x[i]

    return x

def get_covered_positions(df):
    positions = set()
    for i, row in df.iterrows():
        positions.add((row['chromosome'],row['start'],row['end']))
    return positions

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
    for input_file in args.input:
        # read input data from mbtools
        try:
            df = pd.read_csv(input_file, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        # remove positions with zero depth
        df = df[df.num_called_reads > 0]
        
        sample_name.append(get_sample_name(input_file))
        atlas = ReferenceAtlas(args.atlas, get_covered_positions(df))
        t = np.array(df.num_called_reads)
        xhat = np.array(df.modification_frequency)
        m = np.rint((t * xhat))

        # convert to Samples and run
        s = Sample(args.name, xhat, m, t)
        if args.model == 'nnls':
            Y.append(fit_nnls(atlas, s))
        else:
            Y.append(fit_llse(atlas, s, args.epsilon))
    # output
    print("\t".join(['ct'] + sample_name))
    for i, cell_type in enumerate(atlas.get_cell_types()):
        print("\t".join([cell_type] + [str(y[i]) for y in Y]))

if __name__ == "__main__":
    main()
