import pyranges as pr
import pandas as pd
import numpy as np

class ReferenceAtlas:
    """
    Reference atlas class for storing the methylation propensities of each cell type
    """
    def __init__(self, gr):
        """
        Initialize the reference atlas

        :param gr: pyranges object with the following columns: ['Chromosome', 'Start', 'End', 'cell-type1', 'cell-type2', ...]
        :return: None
        """
        self.cpg_ids = [(chrom, start, end) for chrom, start, end in\
                        zip(gr.Chromosome, gr.Start, gr.End)]
        cell_types = list(gr.columns[3:])
        self.K = len(cell_types)
        self.v = {k:list(gr[k]) for k in cell_types}
        self.A = np.array(gr.loc[:, list(cell_types)])

    def get_x(self, sigma):
        """
        Compute the expected methylome by matrix multiplication of the reference atlas and the cell-type proportions

        :param sigma: cell-type proportions
        :return: expected methylome
        """
        return np.dot(self.A, sigma)
        x = np.matmul(self.A, sigma)
        return x

    def get_num_cpgs(self):
        return len(self.cpg_ids)

    def get_cell_types(self):
        return list(self.v.keys())

    def get_num_cell_types(self):
        return len(self.v.keys())

class Sample:
    """
    Sample class for storing the observed methylome
    """
    def __init__(self, name, x_hat, m, t):
        """
        Initialize the sample

        :param name: sample name
        :param x_hat: observed methylome (nnls)
        :param m: observed modified calls (binomial)
        :param t: observed total calls (binomial)
        :return: None
        """
        self.name = name
        self.x_hat = x_hat
        self.m = m
        self.t = t
