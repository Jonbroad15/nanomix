import pyranges as pr
import pandas as pd
import numpy as np
from tools import *

class AtlasMethylome:
    """
    Reference atlas class for storing the methylation propensities of each cell type
    """
    def __init__(self, methylome, atlas, threads=1):
        """
        Take a methylome of read frequencies and turn it into a methylome of region frequencies from regions in the atlas

        Methylome has columns:
        read_name       chromosome      start_position  end_position    alignment_length        strand  mapping_quality total_calls     modified_calls  modification_frequency  cell_type
    e5349968-9ed4-43e4-882b-061ecc7e2683    chr1    1071849 1072048 52806   -       60      16      16      1.00    HCT116
    e5349968-9ed4-43e4-882b-061ecc7e2683    chr1    1066168 1066385 52806   -       60      6       6       1.00    HCT116


        We want to aggregate the columns modified_calls and total_calls based on regions in the first three columns of the atlas file

        :param atlas: Atlas file path
        :param methylome: Methylome file path
        :return: self
        """

        # Rename columns 
        columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                                'start':'Start',
                                'end':'End',
                                'start_position':'Start',
                                'end_position':'End'}

        # Read atlas
        df_atlas = pd.read_csv(atlas, sep='\t').rename(columns=columns)
        df_atlas.drop_duplicates(inplace=True)
        if 'label' in df_atlas.columns: df_atlas.drop('label', axis=1, inplace=True)
        df_atlas.dropna(inplace=True)
        gr_atlas = pr.PyRanges(df_atlas).sort()

        # Read methylomes data from mbtools
        try:
            df = pd.read_csv(methylome, sep='\t').rename(columns=columns)
        except pd.errors.EmptyDataError:
            Exception("Empty methylome file")
        df.dropna(inplace=True)
        gr_sample = pr.PyRanges(df).sort()
        gr = gr_atlas.join(gr_sample, nb_cpu=threads)

        # Check for empty upon join
        cell_types = get_cell_types(atlas)
        if len(gr) == 0:
            Exception("Empty join between atlas and sample. The sample does not overlap with any regions in the atlas.")
        df = gr.df
        df_grouped = df.groupby(['Chromosome', 'Start', 'End'], observed=True)
        df_grouped_sample = df_grouped[['modified_calls', 'total_calls']].sum().sort_index()
        df_grouped_atlas = df_grouped[cell_types].first().sort_index()
        df_join = df_grouped_sample.join(df_grouped_atlas)

        self.t = np.array(df_join.total_calls, dtype=np.float32)
        self.m = np.array(df_join.modified_calls, dtype=np.float32)
        self.x_hat = self.m/self.t

        self.cpg_ids = [(chrom, start, end) for chrom, start, end in df_join.index]
        self.K = len(cell_types)
        self.v = {k:list(df_join[k]) for k in cell_types}
        self.A = np.array(df_join.loc[:, list(cell_types)])

    def get_x(self, sigma):
        """
        Compute the expected methylome by matrix multiplication of the reference atlas and the cell-type proportions

        :param sigma: cell-type proportions
        :return: expected methylome
        """
        return np.dot(self.A, sigma)

    def get_x(self, sigma, p01, p11):
        """
        Compute the expected methylome by matrix multiplication of the reference atlas and the cell-type proportions
        with vectorized p01, p11

        :param sigma: cell-type proportions
        :param p01: p01 vector
        :param p11: p11 vector
        :return: expected methylome
        """
        return np.dot(self.A, sigma*p11) + np.dot(1-self.A, sigma*p01)

    def get_cell_types(self):
        return list(self.v.keys())

    def get_num_cell_types(self):
        return len(self.v.keys())

    def __len__(self):
        return len(self.cpg_ids)

    def __repr__(self):
        return "AtlasMethylome with {} CpGs and {} cell types".format(len(self), self.get_num_cell_types())

