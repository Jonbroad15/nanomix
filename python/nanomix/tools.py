import csv
import numpy as np

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

