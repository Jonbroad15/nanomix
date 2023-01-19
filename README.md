# Nanomix: Methylation Deconvolution of Nanopore Sequencing Data
Methylation deconvolution is the process of determining the proportion of distinct cell types in a complex (hetergeneic) mixture of cell or cell free DNA.
This tool provides suitable models for performing deconvolution on Nanopore sequencing data. In particular our new models account for the non-uniform coverage distribution and high error rate in modified base calling. We also include more typical deconvolution models for deconvolution of bisulfite sequencing data or bead chip arrays.

## Installation
This package is available on PyPI
```
pip install nanomix
```
or alternatively on Conda:
```
conda install nanomix
```

## Usage
To deconvolute a sequencing run, one must simply provide `nanomix` with a methylome. We define a methylome as a `tsv` file with columns `{chr, start, end, total_calls, modified_calls}`. Such a file can be created from a `.bam` file using our associated program, [mbtools](https://github.com/jts/mbtools)
Then the mixture proportion can be found by calling:
```
nanomix deconvolute <METHYLOME>
```

### Reference Atlas
Deconvolution determines the mixture proportion based on methylation propensities of previously resolved sequencing runs of purified reference cells across the genome. This information is collated into an *atlas*. We suggest using the atlas from [Loyfer et. al](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1.full) which we have curated and labelled as `39Bisulfite.tsv` and is set for default. Their tool [wgbstools](https://github.com/nloyfer/wgbs_tools) also provides the means to create an atlas suited to the cell types you are interested in.
```
nanomix deconvolute <METHYLOME> -a <ATLAS>
```

### Model
We provide four deconvolution models

- **llse (default):**   log-likelihood with sequencing errors. Maximize the likelihood of the methylome, atlas and sigma
                    by assuming modification calls follow a binomial distribution. Good for sequencing data with high error
                    rate and non-uniform coverage distribution. (Oxford Nanopore)
- **nnls:**             non-negative least squares. Minimize the squared error between the methylome and what we expect for
                    the methylome (given sigma and the atlas). Recommended for fast deconvolution of methylomes with high
                    coverage. (Methylation Arrays)
- **mmse:**             mixture model with sequencing errors. Also follows a binomial distribution, but softly assigns fragments
                    to cell-types. Optimization uses expectation maximization (slower than above). Recommended for high resolution
                    deconvolution (many cell types) and an atlas with large regions of grouped CpGs.
- **llsp:**             log-likelihood with sequencing perfect. Same as llse, without error modelling. Useful for differentiating the
                    effect of sequencing errors on deconvolution loss and accuracy.
Select a model by:
```
nanomix deconvolute <METHYLOME> -m <MODEL>
```
For more info on other option hparams, run
```
nanomix deconvolute -h
```

### Assign fragments
Our tools also allows you to assign fragments in the methylome to cell types in the atlas based off of the deconvoluted sigma vector.
```
nanomix assign <METHYLOME> -s <SIGMA> 
```
### Simulate 
We provide functionality to simulate methylomes of complex cell mixtures given a `sigma.tsv` file that indicates the cell\_type in the first column and the corresponding proportion in the second column. All the proportions must add up to 1 and the cell-types must be the same as those in the supplied reference atlas. To simulate a methylome:
```
nanomix simulate <SIGMA>
```

### Evaluate
Simulating data provides true cell-type assignments in the last column of the methylome. We can evaluate the performance of a models deconvolution on this methylome. This will output the deconvolution loss (euclidean distance between true and predicted sigma vector) and the read assignment accuracy at confidence levels from 0.5 to 0.9.
```
nanomix evaluate <METHYLOME>
```





