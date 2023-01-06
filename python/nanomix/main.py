
def get_sample_name(s):
    s = s.split('/')[-1]
    # label with coverage level
    # return re.search("\.\d+(\.\d+)*", s)[0][1:]
    return  s.split('.')[0]

def deconvolve(methylomes, atlas, model, p01, p11, random_inits):
    """
    :param methylome: path to tsv file of methylome
    :param atlas: path to tsv file of atlas
    :param model: deconvolution model options: [nnls, llse, llsp, mmse]
    :param p01: sequencing miscall rate
    :param p11: sequencing correct call rate
    :param random_inits: use random initializations for LLSE
    :return: deconvolution results
    """

    if model == 'mmse':
        sigma = fit_mmse(atlas, methylome, p01, p11, stop_thresh=1e-3, max_iter=100, min_proportion=0.01)
    else:
        Y = []
        sample_names = []
        columns={'chromosome':'Chromosome', 'chr':'Chromosome',
                                'start':'Start',
                                'end':'End'}
        df_atlas = pd.read_csv(atlas, sep='\t').rename(columns=columns)
        df_atlas.drop_duplicates(inplace=True)
        if 'label' in df_atlas.columns: df_atlas.drop('label', axis=1, inplace=True)
        df_atlas.dropna(inplace=True)
        gr_atlas = pr.PyRanges(df_atlas).sort()
        for methylome in methylomes:
            # read methylomes data from mbtools
            try:
                df = pd.read_csv(methylome, sep='\t').rename(columns=columns)
            except pd.errors.EmptyDataError:
                continue
            df.dropna(inplace=True)
            gr_sample = pr.PyRanges(df).sort()

            # df_sample = gr_sample.df.groupby(['Chromosome', 'Start', 'End'], as_index=False).sum()
            # gr_sample = pr.PyRanges(df_sample)
            # Init atlas and sample
            gr = gr_atlas.join(gr_sample)
            atlas = ReferenceAtlas(gr.df.loc[:, gr_atlas.columns])
            t = np.array(gr.total_calls, dtype=np.float32)
            m = np.array(gr.modified_calls, dtype=np.float32)


            xhat = m/t
            name = get_sample_name(methylome)
            sample_names.append(name)
            s = Sample(name, xhat, m, t)

            # Run
            if model == 'nnls':
                sigma = fit_nnls(atlas, s)
            elif model == 'llse':
                sigma = fit_llse(atlas, s, p01, p11, random_inits)
            elif model == 'llsp':
                sigma = fit_llsp(atlas, s, p01, p11, random_inits)
            elif model == 'null':
                sigma = fit_uniform(atlas, s)
            else:
                Exception(f"no such model {model}")

            Y.append(sigma)
            # print("name:\t{}".format(name))
            # print("log-likelihood:\t{:.2f}".format(log_likelihood_sequencing_with_errors(atlas, sigma, s, p01, p11)))
            # true_sigma = np.zeros(25)

            # for i, cell_type in enumerate(atlas.get_cell_types()):
                # if cell_type == 'Lung cells':
                    # true_sigma[i] = 0.3
                # elif cell_type == 'Monocytes EPIC':
                    # true_sigma[i] = 0.7
            # print("with true sigma ll:\t{:.6f}".format(log_likelihood_sequencing_with_errors(atlas, true_sigma, s, p01, p11)))

    return Y, sample_names, atlas

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', type=str,
            default='/.mounts/labs/simpsonlab/users/jbroadbent/code/cfdna/nanopore_cfdna/atlases/meth_atlas.csv')
    parser.add_argument('--model', default='llse', type=str, help='deconvolution model options: [nnls, llse, llsp, mmse]')
    parser.add_argument('input', nargs='+',
                        help='reference_modifications.tsv file')
    parser.add_argument('--p01', default=0.05, type=float)
    parser.add_argument('--p11', default=0.95, type=float)
    parser.add_argument('--random_inits', action='store_true')
    args = parser.parse_args()

    Y, sample_names, atlas = deconvolve(args.input, args.atlas, args.model, args.p01, args.p11, args.random_inits)

    if len(Y) < 1: Exception("No output, Deconvolution Failed")
    print("\t".join(['ct'] + sample_names))
    for i, cell_type in enumerate(atlas.get_cell_types()):
        print("\t".join([cell_type] + [str(round(y[i],4)) for y in Y]))



    # print log-likelihood




if __name__ == "__main__":
    main()
