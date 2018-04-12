# archaic_match

Calculate archaic match percent for haplotype windows

## Installation

Due to an [issue with scikit-allel installation](https://github.com/cggh/scikit-allel/issues/177)
numpy must be installed before archaic_match.

```
    pip install numpy
    pip install archaic_match
```

Other requirements
*   sqlite3
*   tabix


## Usage

### Calculate match percent counts to create null distribution database
```
    archaic_match max-match-pct
        --vcf data/simulated_test/null/vcfs/Tenn_nonAfr_*_n1_0.1_n2_0.0.mod.vcf.gz
        --archaic-populations Neand1  # archaic population to match against
        --modern-populations EUR AFR  # populations to generate match percents for
        --chrom-sizes 1000000  # or file with chrom and size columns
        --populations data/simulated_test/null/Tenn.popfile  # file with columns for sample and population
        [--window-size [BP] (default: 50000)]
        [--step-size [BP] (default: 10000)]
        > output/simulated_test/null_tables/afr_eur-vs-neand1/tsv/max_match_pct_counts_all.tsv
```

### Build database from match percent count table(s)
```
    archaic_match build_db
        --match_pct_count output/simulated_test/null_tables/afr_eur-vs-neand1/tsv/max_match_pct_counts_all.tsv
        --db output/simulated_test/null_tables/afr_eur-vs-neand1/max_match_pct_test.db
```

### Calculate match percents and pvalues

Compares each modern haplotype to each archaic haplotype.

For each window (as defined by window_size and step_size), calculates the
frequency of informative sites and the maximum archaic match percent.

Additionally, a empirical pvalue is computed using the database of precomputed
match percents using windows that have an informative site frequency within
`frequency_threshold` of the query window.

```
    archaic_match max-match-pct
        --vcf data/simulated_test/0.1pct/vcfs/Tenn_nonAfr_*_n1_0.1_n2_0.0.mod.vcf.gz
        --archaic-populations Neand1  # archaic population to match against
        --modern-populations EUR AFR  # populations to generate match percents for
        --chrom-sizes=1000000  # or file with chrom and size columns
        --populations data/simulated_test/null/Tenn.popfile  # file with columns for sample and population
        [--window-size [BP] (default: 50000)]
        [--step-size [BP] (default: 10000)]
        --match-pct-database output/simulated_test/null_tables/afr_eur-vs-neand1/max_match_pct.db
        [--frequency-threshold [FLOAT] (default: 0.0001)]
        > output/simulated_test/0.1_pct_pvalues.txt
```
