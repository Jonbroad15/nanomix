import pyranges as pr
import pandas as pd
import numpy as np

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
