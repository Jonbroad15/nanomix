#! /usr/bin/env python
#############################################################
#                                                           #
#  Code adapted from https://github.com/nloyfer/meth_atlas  #
#                                                           #
#############################################################
import numpy as np
import pandas as pd
import os.path as op
import sys
import math
import matplotlib.pylab as plt
import matplotlib.cm
import matplotlib.colors

# Plotting parameters:
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'tab20'         # color map. See https://matplotlib.org/users/colormaps.html
#COLOR_MAP = 'Vega10'

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random

cell_types = ["B-cell", "CD4T-cell", "CD8T-cell", "NK-cell", "T-cell",  "granulocyte",  "breast",  "lung", "lung_alveolar", "monocyte", "neuron", "breast_basal", "breast_luminal", "colon", "hepatocyte", "adipocyte", "bladder", "erythrocyte_progenitor", "bone_osteoblast","kidney", "endothelial_cell",  "lung_endothelial", "kidney_endothelial", "vascular_endothelial", "pancreas_endothelial", "prostate", "epiderminal_keratinocyte", "fallopian", "gallbladder", "gastric","head_neck", "heart_cardiomyocyte", "heart_fibroblasts","left_atrium","oligodendrocyte", "ovary", "pancreas_acinar", "pancreas_alpha", "pancreas_beta", "pancreas_delta", "pancreas_duct","skeletal_muscle", "small_intestine", "smooth_muscle", "thyroid", "upper_GI", "uterus_cervix", "colon_fibroblasts", "dermal_fibroblast", "HSAEC", "HCT116"]

color_map = {}
hatch_map = {}
colors = plt.get_cmap(COLOR_MAP).colors
hatches = [None,'xxx', '///']

for i, cell_type in enumerate(cell_types):
    color_map[cell_type] = colors[i % len(colors)]
    hatch_map[cell_type] = hatches[(i // len(colors)) % len(hatches)]
color_map['other'] = 'black'
hatch_map['other'] = None



def hide_non_blood(df):
    blood_cells = ['tcell', 'erythroblast', 'nkcell', 'bcell', 'progenitor', 'hsc', 'monocyte', 'macrophage', 'eosinophil', 'neutrophil']
    selection = [name not in blood_cells for name in df.index]
    others = df[selection].sum()
    df = df.drop(df.index[tuple([selection])])
    df = df.append(others.rename('other'))
    return df

def aggregate_tissues(df, tissues, name):
    """
    Sum up mixtures of rows in tissues

    :param df: pandas dataframe of mixture proportions
    :tissues: list of tissues to aggregate together
    :return: df with rows aggregated
    """
    tissues = [t for t in tissues if t in df.index]
    # get the rows to aggregate
    rows = [df.loc[tissue] for tissue in tissues]
    # sum them up
    new_row = pd.Series(np.sum(rows, axis=0), name=tissues[0])
    # drop the rows
    df = df.drop(tissues, axis=0)
    # add the new row and give it the name: name
    df = df.append(new_row.rename(name))
    return df


def hide_small_tissues(df, threshold=0.01):
    """
    tissues with very small contribution are grouped to the 'other' category.
    :return: The DataFrame with the new category ('other'),
             where the low-contribution tissues are set to 0.
    """
    others = df[df < threshold].sum()
    df[df < threshold] = 0.0
    # remove rows with all zeros:
    df = df.loc[(df != 0).any(axis=1)]
    df = pd.concat(df, others.rename('other'))
    return df


def gen_bars_colors_hatches(nr_tissues):
    """
    Generate combinations of colors and hatches for the tissues bars
    Every tissue will get a tuple of (color, hatch)
    the last tuple is for the 'other' category, and is always black with no hatch.
    :return: a list of tuples, with length == nr_tissues
    """
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    hatches = [None, 'xxx', '...', 'O', '++'][:nr_tissues // 7]

    nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)
    print(nr_colors, hatches)

    # generate bars colors:
    cmap = matplotlib.cm.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
    colors = [cmap(norm(k)) for k in range(nr_colors)]
    # for c in range(nr_colors): print(norm(c), cmap(norm(c)), cmap(c))

    def get_i_bar_tuple(i):
        color_ind = i % nr_colors
        hatch_ind = int(i // math.ceil(nr_tissues / len(hatches)))
        return colors[color_ind], hatches[hatch_ind]

    colors_hatches_list = [get_i_bar_tuple(i) for i in range(nr_tissues-1)]
    return colors_hatches_list + [((0, 0, 0, 1), None)]

def g(x):
    if x == 'ct':
        return 0
    else:
        return 1/float(x)


def f(x):
    if x == 'ct':
        return 0
    elif 'relapse' in x:
        x_new = "".join(x.split('_')[:2] + ['_'])
        return sum([ord(x0) for x0 in x_new]) + 0.0001
    elif 'BCR' in x:
        return sum([ord(x0) for x0 in x]) + 500
    elif 'add' in x:
        return 1
    else:
        return sum([ord(x0) for x0 in x])

def plot_res(df, outpath, show=False):

    df = hide_small_tissues(df)
    # df = hide_non_blood(df)
    nr_tissues, nr_samples = df.shape

    plt.figure(figsize=FIG_SIZE)
    r = [i for i in range(nr_samples)]
    bottom = np.zeros(nr_samples)
    for i in range(nr_tissues):
        plt.bar(r, list(df.iloc[i, :]),
                edgecolor='white',
                width=0.85,
                label=df.index[i],
                bottom=bottom,
                color=color_map[df.index[i]],
                hatch=hatch_map[df.index[i]])
        bottom += np.array(df.iloc[i, :])

    # Custom x axis
    plt.xticks(r, [str(w)[:NR_CHRS_XTICKS] for w in df.columns], rotation='vertical', fontsize=9)
    plt.xlim(-.6, nr_samples - .4)

    plt.ylabel("Mixture Proportion")
    # Add a legend and a title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title('Deconvolution Results\n'+ op.basename(outpath))

    # adjust layout, save and show
    plt.tight_layout(rect=[0, 0, .83, 1])
    plt.savefig(outpath)
    if show:
        plt.show()

def plot_mixture_proportions(sigma, outpath, group_lung = False):
    """
    Plot the mixture proportions for each sample in a bar chart.
    plots are labeled by their file name

    :param sigma: list of paths to the sigma files. Each containing one deconvolution result.
    :param outpath: path where to save the figure
    """
    # read all input files:
    dfs = []
    for f in sigma:
        df = pd.read_csv(f, sep='\t', index_col=0)
        dfs.append(df)
    # concatenate all dataframes by column, name columns by input file name:
    df = pd.concat(dfs, axis=1, keys=[op.basename(f) for f in sigma])
    # clean up column names
    new_cols = []
    for c in df.columns:
        names = c[0].split('.')
        if names[0] == '0':
            new_cols.append(str(names[1]))
        else:
            new_cols.append(str(names[0]))
    df.columns = new_cols
    # df.columns = [c[0].split('.')[0] for c in df.columns]
    # clean up index names
    df.index = [c.replace(" EPIC", '') for c in df.index]
    # sort columns by name:
    df = df.sort_index(axis=1)
    # Group lung tissues:
    if group_lung:
        df = aggregate_tissues(df, ['lung', 'lung_alveolar', 'lung_endothelial', 'endothelial_cell'], 'lung')
    plot_res(df, outpath)


