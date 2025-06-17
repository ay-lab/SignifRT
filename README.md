# Overview

A Quantitative Framework for the Assessment of RT Changes upon CRISPR-mediated Dissection, from Turner and Hinojosa-Gonzalez et al. (2025). Master transcription factor binding sites constitute the core of early replication control elements. This repository contains all custom functions used in our study to analyze Replication Timing (RT) data.

## Reproducibility

We provide a script that reproduces all main and extended view figures from out study. Simply clone our repository and run the following command after ensuring all required packages are installed in your R environment, and modifying the paths in the example. You can find the necessary data files under `./source_data/` in ([Zenodo](https://zenodo.org/records/15678082)) under DOI 10.5281/zenodo.15678082.

``` bash
ODIR=$(pwd)/output_dir/rep
mkdir -p $ODIR
Rscript  ./reproduce_figs.R -s ./Significance.R -i ./source_data/dppa.tsv -t ./source_data/pre_calc_stats.tsv -o $ODIR
```

## Required R packages

``` r
# Install dependencies
install.packages(c("stringr", "ggplot2", "DescTools", "dplyr", "tidyr","gridExtra","grid")) 
```

## Example Usage: Running on your own RT data.

To employ our AUC based approach to your own data, first ensure that you have installed all dependencies. Then, clone this repository or download the `Significance.R` file directly.

Load your data as a data.frame in R, ensuring that each column corresponds to a different replicate or condition, and that each row corresponds to a genomic bin with the first columns indicating the genomic coordinates. Refer to the toy example below:

``` r
  chr   start   end      c1.2    c2.1    c1.2    c2.2
  <chr> <int> <int>      <dbl>   <dbl>  <dbl>    <dbl>
1 chr1  0     50000       0.51     0.52   0.55     0.54
2 chr1 50000  100000      0.78     0.80   0.82     0.79
3 chr1 100000 150000      0.95     0.95   0.96     0.94
4 chr1 150000 200000      0.81     1.2    0.85     0.88
5 chr1 200000 250000      0.83     0.9    0.87     0.89
```
Note that column names should be named `c.{condition_number}.{replicate_number}`.

You can now run the AUC based analysis using the `Significance()` function:

``` r
source("path/to/Significance.R")
results <- Significance(query='path/to/query.bed',
                        regions=your_data,
                        blacklisted=NA,
                        exp_name='cond1_vs_cond2',
                        type_test='del',
                        size=2000000, 
                        res=50000,
                        n=10000,
                        exclude_query=TRUE,
                        method='area'
)
```
Where the parameters correspond to:

| Parameter | Type | Description |
|:--:|:--:|:-------------------------------------------|
| query | str |Path to a file with genomic coordinates of query regions, or "Dppa" to use the domain used in the paper |
| regions | data.frame | The prepared data.frame with RT values |
| blacklisted | str |Path to a file with genomic coordinates of blacklisted regions that will be excluded from background, or NA if not applicable |
| exp_name | str | Name of the experiment |
| type_test | str |Type of test, either 'del' dor delay, 'adv' for advance or 'chg' if direction does not matter |
| size | int |Size of the genomic bins in bp when building the background model, recommend matching query size |
| res | int |Resolution of the RT data in bp |
| n | int |number of random regions to sample |
| exclude_query | boolean |whether to exclude query regions from the background model |
| method | str |method to use for delta calculation, only 'area' is currently supported |


## About 

