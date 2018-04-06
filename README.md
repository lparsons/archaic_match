# archaic_match

Calculate archaic match percent for haplotype windows

## Installation

Due to an [issue with scikit-allel installation](https://github.com/cggh/scikit-allel/issues/177)
numpy must be installed before archaic_match.

```
    pip install numpy
    pip install archaic_match
```

## Usage


### Calculate maximum archaic match percent for samples

```
    archaic_match max-match-pct
    --vcf-file=[VCF_FILE]
    --archaic-samples=[FILE]
    --modern-samples=[FILE]
    --chrom-sizes=[FILE]
    [--window-size=[BP] (default: 50000)]
    [--step-size=[BP] (default: 10000)]
```

# TODO Ensure informative sites are mixture in modern population

# TODO Check match pct cal`

# TODO Store and lookup haplotype population (Tenn.popfile)

Database for a given archaic population and given modern pop
Lookups done based on specified pop, poplist, or popBySample

For now, using Snakemake to parallelize, could be simpler.


### Generate database tables from text output

Snakemake script for now, but should be command


### Lookup pvalues for samples

Compares each modern haplotype to each archaic haplotype.

For each window (as defined by window_size and step_size), calculates the
frequency of informative sites and the maximum archaic match percent.

Additionally, a emperical pvalue is computed using the database of precomputed
match percents using windows that have an informative site frequency within
`frequency_threshold` of the query window.

```
    archaic_match max-match-pct
        --match-pct-database=[MAX_MATCH_PCT_DATABASE]
        --vcf-file=[VCF_FILE]
        --archaic-samples=[FILE]
        --modern-samples=[FILE]
        --chrom-sizes=[FILE]
        [--window-size=[BP] (default: 50000)]
        [--step-size=[BP] (default: 10000)]
        [--frequency-threshold=[FLOAT] (default: 0.0005)]
```
