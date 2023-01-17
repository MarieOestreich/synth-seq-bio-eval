
<p align="center">
  <img src="https://user-images.githubusercontent.com/50077786/212926291-26e7aea9-cf6e-4464-8f88-440bea2d3a9a.png" width="324" height="192">
</p>

# Biology-based Evaluation of Synthetic Transcriptomics Data
This repo lets you assess if your synthetic transcriptomics dataset is biologically sound. It uses the original train dataset of your generative model as reference and computes the preservation of two common transcriptomics metrics in your synthetic data: Differential expression and gene co-expression.

This repo is part of the **Helmholtz AI project [ProGeneGen](https://progenegen.hmsp.center/)**, a joint project of **[DZNE](https://www.dzne.de/en/)** and **[CISPA](https://cispa.de/en)** to develop generative models for privacy preserving synthetic transcriptomics data.


## How to run it

1. You need your real and synthetic data in the following format: 
    - `.csv` file
    - rows represent samples
    - 1st column contains labels
    - all other columns are genes
2. pass the paths of the two `.csv` files to `PATH_TO_REAL` and `PATH_TO_SYNTH` in the `run.sh`
3. give your run a name that will be used as a file prefix to differentiate plots and reports from different runs (`EXPERIMENT_NAME` in `run.sh`)
4. run `run.sh`

## Outputs

1. a folder named `plots` containing 3 plots per run
    - boxplot of 'false down': expression values of genes the were differentially expressed (down) in the synthetic data, but not in the real data. The boxplot shows their expression in the real vs the synthetic data for all label/class comparisons
    - boxplot of 'false up': expression values of genes the were differentially expressed (up) in the synthetic data, but not in the real data. The boxplot shows their expression in the real vs the synthetic data for all label/class comparisons
    - a plot of the % of maintained co-expressed genes in the synthetic data in comparison to the real data for a variety of correlation cutoffs (based on Pearson Correlation)
2. an html-report summarizing all results including the plots
