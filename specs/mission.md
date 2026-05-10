# Project Mission: Nanomix

## Core Vision
Nanomix is dedicated to advancing the field of epigenetic analysis through high-fidelity methylation deconvolution of Nanopore sequencing data. Our mission is to provide the scientific community with robust, computationally efficient, and mathematically grounded tools to resolve complex cell mixtures from native DNA modifications.

## Target Audience
Nanomix is designed for a diverse community of researchers and developers:
- **Computational Biologists & Bioinformaticians:** Seeking accurate tools for cell-type deconvolution in complex tissue or cfDNA samples.
- **Computer Scientists:** Interested in the application of probabilistic modeling and machine learning to genomic data.
- **Academic & Industry Collaborators:** Looking to integrate state-of-the-art Nanopore methylation analysis into clinical workflows, early cancer screening, or liquid biopsy platforms.

## Technical Innovation
While traditional deconvolution methods (like Non-Negative Least Squares - NNLS) are effective for high-coverage methylation arrays, they often fail to account for the unique characteristics of Nanopore sequencing, such as:
- **Sampling Noise:** Non-uniform coverage across the genome.
- **Asymmetric Error Rates:** Context-specific overcalls and undercalls in modification calling.

Nanomix addresses these challenges by introducing novel likelihood-based algorithms:
- **LLSE (Log-Likelihood with Sequencing Errors):** Directly models Nanopore-specific error rates using binomial success parameters.
- **MMSE (Mixture Model with Sequencing Errors):** Employs Expectation-Maximization (EM) to softly assign individual DNA fragments to cell types.

## Broad Goals
1. **Algorithm Leadership:** Establish new benchmarks for accuracy in low-coverage and high-error sequencing environments.
2. **Accessibility:** Package advanced probabilistic models into an easy-to-use, open-source CLI tool.
3. **Interoperability:** Support diverse reference atlases (array-based and WGBS-based) to bridge the gap between different sequencing technologies.
4. **Foundation for Liquid Biopsy:** Enable the detection of trace amounts of tumor-derived or tissue-specific DNA in cell-free DNA (cfDNA) samples, supporting the next generation of non-invasive diagnostics.
