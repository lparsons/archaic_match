# archaic_match

Calculate archaic match percent for haplotype windows

[![CircleCI](https://circleci.com/gh/lparsons/archaic_match.svg?style=svg)](https://circleci.com/gh/lparsons/archaic_match)

## Installation

Due to an [issue with scikit-allel installation](https://github.com/cggh/scikit-allel/issues/177)
numpy must be installed before archaic_match.

```
    pip install numpy
    pip install git+https://github.com/lparsons/archaic_match
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
        [--overlap-regions [BED_FILE]]
        > output/simulated_test/0.1_pct_pvalues.txt
```

#### Optional region overlap

Specifying a bed file with regions will generate two additional
output columns:

*   `region_overlap_bp`: The number of basepairs of the window that
    overlap an introgressed region.

*   `region_informative_sites`: The number of informative sites in the
    window that overlap an introgressed region.

This is helpful for assessing the amount of a window that in known to be
introgressed in simulated data.

The format of the region bedfile should be a tab delimited file with the
following columns:

1.  Chromsome (sequence id where the region is located)
2.  Start (zero based start position of the region)
3.  End (one-based end position of the region)
4.  Sample \[Optional\] (if specified, only windows associated with this
    sample will be reported as having overlap)

**Note for TreeCalls bedfiles**

The TreeCalls bedfiles output by the *???* simulator should have the haplotype
id converted into the associated sample name (e.g. 220 => msp_110). They
should also be combined into a single bedfile. This can be accomplished using
the included `column_replace` command (e.g.):

```
    column_replace \
        data/simulated_test/0.05_pct/TreeCalls/Tenn_nonAfr_*_n1_0.05_n2_0.0.bed.merged.gz \
        -d data/simulated_test/Tenn_haplotype_to_sample.txt \
        -c 4 \
        | sort -k 1,1 -k 2,2n \
        > data/simulated_test/0/05_pct/TreeCalls/combined_introgressed_regions.bed
```
