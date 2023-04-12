# Nanomix: Methylation Deconvolution of Nanopore Sequencing Data
Methylation deconvolution is the process of determining the proportion of distinct cell types in a complex (hetergeneic) mixture of cell or cell free DNA.
This tool provides suitable models for performing deconvolution on Nanopore sequencing data. In particular our new models account for the non-uniform coverage distribution and high error rate in modified base calling. We also include more typical deconvolution models for deconvolution of bisulfite sequencing data or bead chip arrays.


## Installation
This package is available on PyPI
```
pip install nanomix
```
Alternatively you can install from source. Installing from source requires [maturin](https://github.com/PyO3/maturin)
```
pip install maturin
git clone https://github.com/Jonbroad15/nanomix.git
cd nanomix
maturin develop
```

## Usage
### Reference Atlas
Deconvolution determines the mixture proportion based on methylation propensities of previously resolved sequencing runs of purified reference cells across the genome. This information is collated into an *atlas*. We suggest using the atlas from [Loyfer et. al](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1.full) which we have curated and labelled as `39Bisulfite.tsv` and is set for default. Their tool [wgbstools](https://github.com/nloyfer/wgbs_tools) also provides the means to create an atlas suited to the cell types you are interested in.

### Deconvolution
To deconvolute a sequencing run, one must simply provide `nanomix` with a methylome. We define a methylome as a `tsv` file with columns `{chr, start, end, total_calls, modified_calls}`. Such a file can be created from a `.bam` file using our associated program, [mbtools](https://github.com/jts/mbtools)
```
mbtools region-frequency -r ATLAS.tsv SAMPLE.bam > METHYLOME.tsv
```
Then the mixture proportion can be found by calling:
```
nanomix deconvolute -a ATLAS.tsv METHYLOME.tsv
```

### Model
We provide four deconvolution models

- **llse (default):**   log-likelihood with sequencing errors. Maximize the likelihood of sigma
                    by assuming modification calls follow a binomial distribution. Good for sequencing data with high error
                    rate. (Oxford Nanopore)
- **nnls:**             non-negative least squares. Minimize the squared error between the methylome and what we expect for
                    the methylome (given sigma and the atlas). Recommended for fast deconvolution of methylomes with high
                    coverage. (Methylation Arrays)
- **llsp:**             log-likelihood with sequencing perfect. Same as llse, without error modelling. Useful for differentiating the
                    effect of sequencing errors on deconvolution loss and accuracy.
- **mmse:**             mixture model with sequencing errors. Also follows a binomial distribution, but softly assigns fragments
                    to cell-types. Optimization uses expectation maximization (slower than above). Recommended for high resolution
                    deconvolution (many cell types) and an atlas with large regions of grouped CpGs.
Select a model by:
```
nanomix deconvolute -m MODEL METHYLOME.tsv 
```
The **mmse model is distinct** in that it works by assigning reads to cell types. To this effect, one would need a methylome where every row represents a read and columns contain `{read_id, chr, start, end, total_calls, modified_calls}`, this also be constructed from a `.bam` file with [mbtools](https://github.com/jts/mbtools)
```
mbtools read-frequency SAMPLE.bam > METHYLOME.tsv
nanomix deconvolute -m mmse METHYLOME.tsv
```
For more info on other option hparams, run
```
nanomix deconvolute -h
```

### Assign fragments
Our tools also allows you to assign fragments in the methylome to cell types in the atlas based off of the deconvoluted sigma vector.
```
nanomix assign -s SIGMA.tsv METHYLOME.tsv 
```
### Simulate 
We provide functionality to simulate methylomes of complex cell mixtures given a `sigma.tsv` file that indicates the cell\_type in the first column and the corresponding proportion in the second column. All the proportions must add up to 1 and the cell-types must be the same as those in the supplied reference atlas. To simulate a methylome:
```
nanomix simulate -a ATLAS.tsv SIGMA.tsv
```

### Evaluate
Simulating data provides true cell-type assignments in the last column of the methylome. We can evaluate the performance of a models deconvolution on this methylome. This will output the deconvolution loss (euclidean distance between true and predicted sigma vector) and the read assignment accuracy at confidence levels from 0.5 to 0.9.
```
nanomix evaluate -a ATLAS.tsv METHYLOME.tsv
```

### Plot
You can plot a list of deconvolution mixtures by providing them to the plot function. This will produce a stacked bar plot.
```
nanomix plot -o NAME.png *sigma.tsv
```
![exampledeconvplot](Images/example_deconvolution_plot.png)




