import numpy as np
import csv
from collections import Counter
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

from _nanomix import *
from models import fit_model, fit_mmse, log_likelihood_sequencing_with_errors, log_likelihood_sequencing_perfect
from atlas import AtlasMethylome
from tools import *
from plot import *

def plot(sigma, outpath, chart, group_lung=False, cell_types=[]):
    """
    Wrapper function to choose plotting function

    :param sigma:
    :param outpath:
    """
    if chart == 'bar':
        plot_mixture_proportions(sigma, outpath, group_lung=group_lung)
    elif chart == 'scatter':
        plot_scatter(sigma, outpath, cell_types=cell_types)

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
    sigma = get_sigma_init(sigma, get_cell_types(atlas))
    generate_methylome(atlas, sigma, coverage, region_size, p01, p11)

def evaluate(methylome, atlas, sigma,
                p01=1.,
                p11=0.,
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


    # Convert counts to proportions
    total_fragments = sum(cell_type_fragment_counts.values())
    true_cell_type_proportions = {cell_type : count / total_fragments for cell_type, count in cell_type_fragment_counts.items()}
    true_sigma = np.array([true_cell_type_proportions[cell_type] for cell_type in cell_types])

    # open sigma.tsv file and load it
    sigma = get_sigma_init(sigma, get_cell_types(atlas))
    deconvolution_loss = np.linalg.norm(true_sigma - sigma)
    print(f"Deconvolution loss:\t{deconvolution_loss}")

    # Assign cell types
    model = MMSE(methylome, atlas, sigma, p01=p01, p11=p11, concentration=1/len(cell_types))
    cell_types += ['unassigned']
    for threshold in [.5, .6, .7, .8, .9]:
        enumerated_assignments = model.assign_fragments_t(threshold)
        assignments = [cell_types[assignment] for assignment in enumerated_assignments]
        correct = 0
        for true_assignment, assignment in zip(true_assignments, assignments):
            if 'PBMC' in true_assignment and assignment in 'monocyte erythrocyte_progenitor T-cell B-cell NK-cell CD4T-cell CD8T-cell':
                correct += 1
            elif 'CliveOME' in true_assignment and assignment in 'monocyte erythrocyte_progenitor T-cell B-cell NK-cell granulocyte CD4T-cell CD8T-cell':
                correct += 1
            elif true_assignment == assignment:
                correct += 1
        accuracy = correct / len(assignments)
        print(f"Accuracy at confidence {threshold}:\t{accuracy}")

def transform_atlas(atlas, p01, p11, threads=1):
    """
    Transform an atlas to include error parameters fixed in.
    Here p01, p11 are defined for every region in the atlas and every cell_type
    such that atlas, p01, p11 are paths to tsv files with the same column and rows

    :param atlas: path to tsv of atlas file
    :param p01: path to tsv of p01 file
    :param p11: path to tsv of p11 file
    """
    cell_types = get_cell_types(atlas)
    # columns
    columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                            'start':'Start',
                            'end':'End',
                            'start_position':'Start',
                            'end_position':'End'}
    # load atlas, p01, p11 into pyranges
    atlas_pr = pd.read_csv(atlas, sep='\t').rename(columns=columns)
    p01_pr = pr.PyRanges(pd.read_csv(p01, sep='\t').rename(columns=columns))
    p11_pr = pr.PyRanges(pd.read_csv(p11, sep='\t').rename(columns=columns))

    # join p01, p11 to atlas
    pr = atlas_pr\
        .join(p01_pr, suffix=_p01, nb_cpu=threads)\
        .join(p11_pr, suffix=_p11, nb_cpu=threads)\
        .sort(nb_cpu=threads)

    # add columns for each cell type if they do not exist
    # i.e. if lung_p01 not in columns duplicate default_p01 and call it lung_p01
    for cell_type in cell_types:
        if cell_type + '_p01' not in pr.columns:
            pr[cell_type + '_p01'] = pr['default_p01']
        if cell_type + '_p11' not in pr.columns:
            pr[cell_type + '_p11'] = pr['default_p11']

    # make a numpy matrix of the original elements from atlas_pr
    A = pr.as_df()[cell_types].to_numpy()
    P01 = pr.as_df()[[f'{cell_type}{_p01}' for cell_type in cell_types]].to_numpy()
    P11 = pr.as_df()[[f'{cell_type}{_p11}' for cell_type in cell_types]].to_numpy()

    # Do element wise multiplication of A and P11
    # Do element wise multiplication of (1-A) and P01
    # Add the two matrices together
    transformed_atlas = A * P11 + (1 - A) * P01

    # write to stdout with locus and cell_type columns
    print('chromosome\tstart\tend\t' + '\t'.join(cell_types))
    for i, row in enumerate(transformed_atlas):
        row_as_string = '\t'.join([str(x) for x in row])
        print(f"{pr.as_df().iloc[i]['Chromosome']}\t{pr.as_df().iloc[i]['Start']}\t{pr.as_df().iloc[i]['End']}\t{row_as_string}")



def log_likelihood(methylome, atlas, sigma, p01, p11, model):
    """
    Return the log likelihood for a given solution

    :param methylome: Path to tsv file of methylome
    :param atlas: Path to tsv file of atlas
    :param model: Deconvolution model options: [nnls, llse, llsp, mmse]
    :param p01: Sequencing miscall rate
    :param p11: Sequencing correct call rate
    :param sigma: path to sigma file
    :return: None
    """
    atlas_path = atlas
    sigma = get_sigma_init(sigma, get_cell_types(atlas))
    ll = 'none'
    if model == 'llse':
        data = AtlasMethylome(atlas_path, methylome)
        # print("Total number of methylated calls: {}".format(np.sum(s.m)))
        # print("Total number of unmethylated calls: {}".format(np.sum(s.t-s.m)))
        # print("atlas shape: {}".format(atlas.A.shape))
        ll = log_likelihood_sequencing_with_errors(data, sigma, p01, p11)
    elif model == 'mmse':
        mmse = MMSE(methylome, atlas_path, sigma, p01, p11, 0.1)
        ll = mmse.log_likelihood()
    elif model == 'llsp':
        data = AtlasMethylome(atlas_path, methylome)
        ll = log_likelihood_sequencing_perfect(data, sigma, p01, p11)
    print("name\tvalue")
    print("Atlas\t{}".format(atlas_path.split('/')[-1].replace('.tsv', '')))
    print("Model\t{}".format(model))
    print("p01\t{}".format(p01))
    print("p11\t{}".format(p11))
    print("Methylome\t{}".format(methylome.split('/')[-1].replace('.tsv', '')))
    print("Log-likelihood\t{}".format(ll))


def deconvolute(methylome, atlas, model,
                p01=0.,
                p11=1.,
                n_trials=10,
                threads=1,
                concentration=None,
                nnls_init=False,
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

    # if p01 is a string, then it is filepath
    # first try and parse it into a float and then if it doesnt work then it is a filepath
    cell_types = get_cell_types(atlas)
    p01 = get_vectorized_error_param(p01, cell_types)
    p11 = get_vectorized_error_param(p11, cell_types)



    if concentration is None:
        concentration = 1/len(get_cell_types(atlas))
    # Run
    if model == 'mmse':
        if nnls_init:
            sigma = fit_model(methylome, atlas, 'nnls', p01, p11, concentration=concentration, threads=threads)
        else:
            sigma = get_sigma_init(sigma_init, get_cell_types(atlas), concentration=concentration)
        cell_type_proportions = fit_mmse(methylome, atlas, sigma, p01, p11, stop_thresh, max_iter, min_proportion, concentration)
    else:
        sigma = fit_model(methylome, atlas, model, p01, p11,
                          n_trials=n_trials,
                          threads=threads,
                          concentration=concentration,
                          nnls_init=nnls_init)
        if min_proportion > 0.0:
            sigma[sigma < min_proportion] = 0.0
            sigma /= np.sum(sigma)
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


