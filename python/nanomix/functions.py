import numpy as np
import csv
from collections import Counter
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

from _nanomix import *
from models import fit_model, fit_mmse
from atlas import ReferenceAtlas, Sample
from tools import *

def simulate(atlas, sigma, coverage, region_size, p01, p11):
    """
    Simulate a methylome from a reference atlas at the given cell type proportion
    Wrapper function to call Rust code from Python

    :param atlas: path to reference atlas
    :param sigma: path to tsv file containing cell type proportions
    :param coverage: coverage of simulated reads
    :param region_size: number of CpGs in each region
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :return: None
    """
    generate_methylome(atlas, sigma, coverage, region_size, p01, p11)

def evaluate(methylome, atlas, model,
                p01=0.,
                p11=1.,
                sigma_init='null',
                max_iter=10,
                min_proportion=0.01,
                stop_thresh=1e-3):
    """
    Evaluate the performance of a model on a simulated methylome
    Wrapper function to call Rust code from Python

    :param methylome: path to simulated methylome
    :param atlas: path to reference atlas
    :param model: model to evaluate
    :param p01: nanopore miscall rate
    :param p11: nanopore correct call rate
    :param sigma_init: model to initialize sigma (mmse only)
    :param max_iter: maximum number of iterations (mmse only)
    :param min_proportion: minimum proportion of a cell type (mmse only)
    :param stop_thresh: stopping threshold (mmse only)
    :return: None
    """
    # Determine cell_type proportions from methylome
    cell_types = get_cell_types(atlas)
    cell_type_fragment_counts = Counter(cell_types)
    true_assignments = []

    # Open methylome and count 
    with open(methylome, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        cell_type_idx = header.index('cell_type')
        for row in reader:
            cell_type_fragment_counts[row[cell_type_idx]] += 1
            true_assignments.append(row[cell_type_idx])

    # Enumerate true_assignments
    enumerated_true_assignments = [cell_types.index(cell_type) for cell_type in true_assignments]

    # Convert counts to proportions
    total_fragments = sum(cell_type_fragment_counts.values())
    true_cell_type_proportions = {cell_type : count / total_fragments for cell_type, count in cell_type_fragment_counts.items()}
    true_sigma = np.array([true_cell_type_proportions[cell_type] for cell_type in cell_types])

    # Run deconvolution
    if model == 'mmse':
        sigma = fit_model(methylome, atlas, sigma_init, p01, p11)
        cell_type_proportions = fit_mmse(atlas, methylome, sigma, p01, p11, stop_thresh, max_iter, min_proportion,
                                         true_sigma=true_sigma, true_assignments=enumerated_true_assignments)
    else:
        cell_type_proportions = deconvolute(methylome, atlas, model, p01, p11, sigma_init, max_iter, min_proportion, stop_thresh, print_output=False)

    # Compute deconvolution loss
    sigma = np.array([cell_type_proportions[cell_type] for cell_type in cell_types])
    deconvolution_loss = np.linalg.norm(true_sigma - sigma)
    print(f"Deconvolution loss: {deconvolution_loss}")

    # Assign cell types
    cell_types += ['unassigned']
    model = MMSE(methylome, atlas, sigma, p01=0., p11=1.)
    for threshold in [.5, .6, .7, .8, .9]:
        enumerated_assignments = model.assign_fragments_t(threshold)
        assignments = [cell_types[assignment] for assignment in enumerated_assignments]
        accuracy = sum([1 if true_assignment == assignment else 0 for true_assignment, assignment in zip(true_assignments, assignments)]) / len(assignments)
        print(f"Accuracy at confidence {threshold}: {accuracy}")

def deconvolute(methylome, atlas, model,
                p01=0.,
                p11=1.,
                sigma_init='null',
                max_iter=10,
                min_proportion=0.01,
                stop_thresh=1e-3,
                print_output=True):
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


    # Run
    if model == 'mmse':
        sigma = fit_model(methylome, atlas, sigma_init, p01, p11)
        cell_type_proportions = fit_mmse(atlas, methylome, sigma, p01, p11, stop_thresh, max_iter, min_proportion)
    else:
        sigma = fit_model(methylome, atlas, model, p01, p11)
        cell_type_proportions = {cell_type: proportion for cell_type, proportion in zip(get_cell_types(atlas), sigma)}

    # return deconvolution results
    if print_output:
        print("cell_type\tproportion")
        for cell_type, proportion in cell_type_proportions.items():
            print(f"{cell_type}\t{proportion}")

    return cell_type_proportions

def assign_fragments(methylome, atlas, sigma,
                     p01=0., p11=1.,
                     threshold=0.5,
                     print_output=True):
    """
    Assign fragments in the methylome to cell types
    Initialize MMSE model to use maximum responsibility (gamma) assignments

    :param atlas: path to reference atlas
    :param methylome: path to methylome
    :param sigma: path to tsv file containing cell type proportions
    :param threshold: confidence threshold for assigning fragments to cell types
    :param print_output: print output to stdout
    :return: assignments vector in order of fragments in the methylome
    """
    # Need to ensure that the cell types are in the same order as the atlas
    cell_types = get_cell_types(atlas)

    # get sigma vector
    with open(sigma, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        cell_type_proportions = {row[0] : float(row[1]) for row in reader}

    # get sigma vector in the same order as the atlas
    sigma = np.array([cell_type_proportions[cell_type] for cell_type in cell_types])

    model = MMSE(methylome, atlas, sigma, p01=0., p11=1.)
    enumerated_assignments = model.assign_fragments_t(threshold)

    # map enumerated assignments to a vector of (string) cell type assignments
    cell_types += ['unassigned']
    cell_type_assignments = [cell_types[i] for i in enumerated_assignments]

    # return cell type assignments to standard output
    if print_output:
        with open(methylome, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            print('\t'.join(header + ['cell_type']))
            for row, cell_type in zip(reader, cell_type_assignments):
                print('\t'.join(row + [cell_type]))

    return cell_type_assignments


