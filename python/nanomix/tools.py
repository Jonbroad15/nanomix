import csv
import numpy as np
import pandas as pd
import pyranges as pr


def get_vectorized_error_param(file, cell_types):
    """
    # we expect a tsv file with first column cell_type and the second column p01
    # open the file and read p01 values into a vector, ordered by cell_types in the atlas
    # order them by cell_types in the atlas
    """
    try:
        error_param = float(file)
        return np.array([error_param]*len(cell_types))
    except ValueError:
        with open(file, 'r') as f:
            # read tsv into a dictionary mapping cell_type to error_param
            error_param_dict = {row[0]: float(row[1]) for row in csv.reader(f, delimiter='\t')}
        error_param = []
        for cell_type in cell_types:
            if cell_type in error_param_dict:
                error_param.append(error_param_dict[cell_type])
            else:
                error_param.append(error_param_dict['default'])
        return np.array(error_param)

def eq_constraint(x):
    return 1 - np.sum(x)

def get_cell_types(atlas):
    """
    Open atlas file and read the header to get the cell types order

    :param atlas: path to reference atlas
    :return: vector of cell types in same order as the atlas
    """
    with open(atlas, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        cell_types = header[3:]
    return cell_types

def get_sigma_init(sigma, cell_types, concentration=1.):
    """
    Get sigma vector from sigma_init file in the same order as the atlas

    :param sigma: path to sigma_init file
    :param cell_types: vector of cell types in same order as the atlas
    :return: vector of cell-type proportions
    """
    if sigma == 'null':
        return np.array([concentration] * len(cell_types))
    else:
        with open(sigma, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            cell_type_proportions = {row[0]: float(row[1]) for row in reader}
        sigma_init = np.array([cell_type_proportions[cell_type] for cell_type in cell_types])
        return sigma_init









