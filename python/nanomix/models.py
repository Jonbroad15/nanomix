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
from multiprocessing import Queue
import threading
from functools import partial

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

from _nanomix import *
from atlas import AtlasMethylome
from tools import *

def fit_model(methylome, atlas, model, p01, p11,
              threads=1, n_trials=5,
              nnls_init=False, concentration=1.0):
    """
    Wrapper function to select model for deconvolution
    mmse model is initalized seperately with Rust bindings
    methylome and data are joined on covered regions before deconvolution

    :param methylome: path to simulated methylome
    :param data_path: path to reference data
    :param model: model to fit
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: cell-type proportions
    """
    data = AtlasMethylome(methylome, atlas, threads=threads)
    if model == 'nnls':
        sigma = fit_nnls(data)
    elif model == 'llse':
        sigma = fit_llse_parallel(data, p01, p11, n_trials, threads, concentration, nnls_init)
    elif model == 'llsp':
        sigma = fit_llsp(data)
    elif model == 'null':
        sigma = fit_uniform(data)
    else:
        raise ValueError(f"no such model: {model}. Choose from [nnls, llse, llsp, mmse]")
    return sigma

def log_likelihood_sequencing_perfect(data, sigma):
    """
    Compute the log-likelihood of the llsp model

    :param data: reference data
    :param sigma: cell-type proportions
    :param data: data methylome
    :return: log-likelihood
    :rtype: float
    """
    sigma_t = sigma.reshape( (data.K, 1) )
    x = np.clip(np.ravel(data.get_x(sigma_t)), 0, 1.0)
    b =  binom.logpmf(data.m, data.t, x)
    return np.sum(b)

def fit_llsp(data, n_trials=10):
    """
    fit the log-likelihood sequencing perfect model

    :param data: reference data
    :param data: data methylome
    :param n_trials: number of trials to run with random initalizations of sigma
    :return: cell-type proportions
    :rtype: np.array
    """
    f = lambda x: -1 * log_likelihood_sequencing_perfect(data, x)
    bnds = [ (0.0, 1.0) ] * data.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / data.K ] * data.K)
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

def log_likelihood_sequencing_with_errors(data, sigma, p01, p11):
    """
    Compute the log-likelihood of the llse model

    :param data: reference data
    :param sigma: cell-type proportions
    :param data: data methylome
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: log-likelihood
    :rtype: float
    """

    sigma_t = sigma.reshape( (data.K, 1) )

    # if p01 is a vector then reshape if not then make it a vector
    if isinstance(p01, np.ndarray):
        p01 = p01.reshape( (data.K, 1) )
    else:
        p01 = np.array([p01] * data.K).reshape( (data.K, 1) )
    if isinstance(p11, np.ndarray):
        p11 = p11.reshape( (data.K, 1) )
    else:
        p11 = np.array([p11] * data.K).reshape( (data.K, 1) )
    p = np.clip(np.ravel(data.get_x(sigma_t, p01, p11)), 0, 1.0)
    # breakpoint()
    # x = np.clip(np.ravel(np.dot(data.A, sigma_t)), 0, 1.0)
    # p = x * p11 + (1-x) * p01
    b =  binom.logpmf(data.m, data.t, p)
    binomial_coef = sum([math.log(math.comb(int(t), int(m))) for m,t in zip(data.m, data.t)])

    return np.sum(b) #- binomial_coef

def fit_uniform(K):
    """
    fit the uniform model (null model)

    :param K: number of cell-types
    :return: cell-type proportions
    :rtype: np.array
    """
    return np.array([1.0 / K ] * K)

class my_thread(threading.Thread):
    """
    thread class for parallel processing of fit_llse
    """
    def __init__(self, data, p01, p11, queue):
        super(my_thread, self).__init__()
        self.data = data
        self.p01 = p01
        self.p11 = p11
        self.queue = queue
        self.res = []
    def run(self):
        alpha = self.queue.get()
        f = lambda x: -1 * log_likelihood_sequencing_with_errors(self.data, x, self.p01, self.p11)
        bnds = [ (0.0, 1.0) ] * self.data.K
        cons = ({'type': 'eq', 'fun': eq_constraint})
        self.res.append(minimize(f, alpha, method='SLSQP', options={'maxiter': 200, 'disp':False}, bounds=bnds, constraints=cons))


def fit_llse_parallel(data, p01, p11, n_trials, threads, concentration, init_nnls):
    alpha = np.array([concentration] * data.K)
    initializations = dirichlet.rvs(alpha, size=n_trials).tolist()
    if init_nnls:
        initializations.append(fit_nnls(data))

    queues = [Queue() for i in range(threads)]
    res = []
    for i in range(n_trials):
        queues[i % threads].put(initializations[i])
    threads = [my_thread(data, p01, p11, queues[i]) for i in range(threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
        res += t.res

    best_ll = np.inf
    best_sol = None
    for r in res:
        ll = r.get("fun")
        if ll < best_ll:
            best_ll = ll
            best_sol = r
    return best_sol.x/np.sum(best_sol.x)

def fit_llse(data, p01, p11):
    """
    fit the log-likelihood sequencing with errors model

    :param data: reference data
    :param data: data methylome
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :param n_trials: number of trials to run with random initalizations of sigma
    :return: cell-type proportions
    :rtype: np.array
    """
    f = lambda x: -1 * log_likelihood_sequencing_with_errors(data, x, p01, p11)
    bnds = [ (0.0, 1.0) ] * data.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    alpha = np.array([ 1.0 / data.K ] * data.K)
    res = minimize(f, alpha, method='SLSQP', options={'maxiter': 100, 'disp':False}, bounds=bnds, constraints=cons)
    return res.x/np.sum(res.x)

def fit_nnls(data):
    """
    fit the non-negative least squares model

    :param data: reference data
    :param data: data methylome
    :return: cell-type proportions
    :rtype: np.array
    """

    # add sum=1 constraint
    t = np.array([1.0] * data.K).reshape( (1, data.K) )
    A = np.append(data.A, t, axis=0)
    b = np.append(data.x_hat, [1.0], axis=0)
    res = nnls(A, b)
    return res[0]/np.sum(res[0])

def fit_mmse(methylome, atlas, sigma, p01, p11, stop_thresh, max_iter, min_proportion, concentration,
             true_sigma=None, true_assignments=None):
    """
    Fit mixture model with sequencing errors (MMSE) onto data with reference data.
    Wrapper function to access Rust implementation.

    :param methylome: path to sample methylome to fit
    :param atlas: path to reference atlas
    :param p01: sequencing miscall rate
    :param p11: sequencing correct call rate
    :param stop_thresh: threshold for stopping iterations
    :param max_iter: maximum number of iterations
    :param min_proportion: minimum proportion of a cell type in the mixture
    :return: fitted mixture model
    """
    # initialize MMSE model
    mmse = MMSE(methylome, atlas, sigma, p01, p11, concentration)
    if true_sigma is not None and true_assignments is not None:
        mmse.evaluate(stop_thresh, max_iter, min_proportion, true_sigma, true_assignments)
    else:
        mmse.optimize(stop_thresh, max_iter, min_proportion)
    return mmse.cell_type_proportions()

