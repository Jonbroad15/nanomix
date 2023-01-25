import csv
import numpy as np
import pandas as pd
import pyranges as pr
from atlas import ReferenceAtlas, Sample

def load_atlas(atlas_path, methylome):
    """
    Load atlas and methylome data. Perform join on the probes.

    :param atlas: Atlas file path
    :param methylome: Methylome file path
    :return: (atlas, sample) tuple
    """
    # load atlas
    columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                            'start':'Start',
                            'end':'End'}
    df_atlas = pd.read_csv(atlas_path, sep='\t').rename(columns=columns)
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

    # Join atlas and sample
    gr = gr_atlas.join(gr_sample)
    # Check for empty upon join
    if len(gr) == 0:
        Exception("Empty join between atlas and sample. The sample does not overlap with any regions in the atlas.")
    atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])
    t = np.array(gr.total_calls, dtype=np.float32)
    m = np.array(gr.modified_calls, dtype=np.float32)
    xhat = m/t
    sample = Sample('methylome', xhat, m, t)

    return atlas, sample

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

