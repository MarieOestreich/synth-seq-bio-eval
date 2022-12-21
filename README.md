# synth-seq-bio-eval
Biology-based evaluation of synthetic transcriptomic data


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