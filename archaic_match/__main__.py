#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" archaic_match

Calculate archaic match percent for haplotype windows
"""

import argparse
import logging
import allel
import numpy
import subprocess
import sys
import functools
import glob
from .classmodule import Window
from .funcmodule import get_samplename_list
from .funcmodule import get_chrom_sizes
from .funcmodule import generate_windows
from .funcmodule import get_sample_populations


__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2018, Lance Parsons"
__license__ = "MIT https://opensource.org/licenses/MIT"


def main():

    # Top level parser
    parent_parser = argparse.ArgumentParser(
        description="Calculate archaic match percent for haplotype windows",
        epilog="As an alternative to the commandline, params can be placed in "
        "a file, one per line, and specified on the commandline like "
        "'%(prog)s @params.conf'.",
        fromfile_prefix_chars='@',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True)

    # Input file options
    input_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    vcf_group = input_parser.add_argument_group('VCF Files')
    vcf_group.add_argument("--vcf",
                           required=True,
                           nargs="+",
                           help="VCF file to parse",
                           metavar="VCF_FILE")

    population_group = input_parser.add_argument_group('Population options')
    population_group.add_argument("--populations",
                                  required=True,
                                  help="Tab delimited file with columns "
                                  "sample, population, and super-population",
                                  metavar="FILE")
    population_group.add_argument("--archaic-populations",
                                  required=True,
                                  nargs='+',
                                  help="List of archaic populations",
                                  metavar="FILE")
    population_group.add_argument("--modern-populations",
                                  required=True,
                                  nargs='+',
                                  help="List of modern populations",
                                  metavar="FILE")
    population_group.add_argument("--chrom-sizes",
                                  required=True,
                                  help="Tab delimited file with seqid\tlength "
                                  "or INT specifying default chromsome length",
                                  metavar="FILE")

    pvalue_group = input_parser.add_argument_group('Pvalue computation')
    pvalue_group.add_argument("--match-pct-database",
                              required=False,
                              help="SQLite database containting distribution "
                              "of archaic match percents for various windows "
                              "from population with no intregression",
                              metavar="SQLITE_FILE")
    pvalue_group.add_argument("--frequency-threshold",
                              required=False,
                              help="Use windows from the database that have "
                              "informative site frequencies within THRESHOLD "
                              "of the query windows informative site "
                              "frequency to determine the null distribution",
                              default=0.0001,
                              type=float,
                              metavar="THRESHOLD")

    # Window size options
    window_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    window_group = window_parser.add_argument_group('Window size options')
    window_group.add_argument("--window-size",
                              help=("Size of windows to calculate match "
                                    "percent"),
                              default="50000",
                              type=int)
    window_group.add_argument("--step-size",
                              help=("Number of basepairs to step to calculate "
                                    "sliding window"),
                              default="10000",
                              type=int)

    # Output file options
    # output = parser.add_argument_group('Output files')
    # output.add_argument("windows_file",
    #                     help="Window information (informative site "
    #                     "frequency)",
    #                     metavar="windows.bed")

    # Other options
    common_options = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    common_options.add_argument("-v", "--verbose",
                                help="increase output verbosity",
                                action="store_true")
    common_options.add_argument("--version", action="version",
                                version="%(prog)s " + __version__)

    subparsers = parent_parser.add_subparsers(
        help="commands")
    parser_max_match_pct = subparsers.add_parser(
        'max-match-pct',
        parents=[input_parser, window_parser, common_options])
    parser_max_match_pct.set_defaults(func=max_match_pct)
    # parser_max_match_pct_pvalue = subparsers.add_parser(
    #     'max-match-pct-pvalue',
    #     parents=[input_parser, window_parser])
    # parser_max_match_pct_pvalue.set_defaults(func=max_match_pct_pvalue)

    # Parse arguments
    args = parent_parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    if args.func:
        args.func(args)


@functools.lru_cache(maxsize=None)
def max_match_pct_pvalue(informative_site_frequency, population,
                         max_match_pct, dbconn, threshold):
    '''Calculate emperical pvalue from windows in database'''
    lower_threshold = informative_site_frequency - threshold
    upper_threshold = informative_site_frequency + threshold
    n = windows_within_isf_threshold(
        informative_site_frequency=informative_site_frequency,
        population=population,
        threshold=threshold,
        dbconn=dbconn)
    c = dbconn.cursor()
    c.execute('''SELECT count(*) from MAX_MATCH_PCT
              where informative_site_frequency >= ?
              and informative_site_frequency <= ?
              and population = ?
              and max_match_pct >= ?''',
              (lower_threshold, upper_threshold, population, max_match_pct))
    s = float(c.fetchone()[0])
    pvalue = (s + 1) / (n + 1)
    return pvalue


@functools.lru_cache(maxsize=None)
def windows_within_isf_threshold(informative_site_frequency, population,
                                 threshold, dbconn):
    '''Return number of windows for a given population within the informative
    site frequency threshold'''
    lower_threshold = informative_site_frequency - threshold
    upper_threshold = informative_site_frequency + threshold
    c = dbconn.cursor()
    c.execute('''SELECT count(*) from MAX_MATCH_PCT
              where informative_site_frequency >= ?
              and informative_site_frequency <= ?
              and population = ?''',
              (lower_threshold, upper_threshold, population))
    n = float(c.fetchone()[0])
    return n


def max_match_pct(args):

    sample_populations = get_sample_populations(args.populations)
    logging.debug("Sample populations: {}".format(sample_populations))

    archaic_sample_list = get_samplename_list(args.archaic_populations,
                                              sample_populations)
    logging.debug("Archaic sample list: {}".format(archaic_sample_list))

    modern_sample_list = get_samplename_list(args.modern_populations,
                                             sample_populations)
    logging.debug("Modern sample list: {}".format(modern_sample_list))

    if args.chrom_sizes.isdigit():
        default_chromsize = int(args.chrom_sizes)
        chrom_sizes = dict()
    else:
        default_chromsize = None
        chrom_sizes = get_chrom_sizes(args.chrom_sizes)
    logging.debug("Chromosome sizes: {}".format(chrom_sizes))

    if args.match_pct_database:
        import sqlite3
        dbconn = sqlite3.connect(args.match_pct_database)

    vcf_filelist = list()
    for vcf_file_pattern in args.vcf:
        vcf_filelist.extend(glob.glob(vcf_file_pattern))
    logging.debug("VCF File list:\n{}".format(vcf_filelist))
    for vcf_file in vcf_filelist:
        chromsomes_in_vcf = subprocess.run(['tabix', '-l', vcf_file],
                                           stdout=subprocess.PIPE,
                                           universal_newlines=True)
        # logging.debug(chromsomes_in_vcf.stdout)
    #    with open(args.windows_file, 'w') as windows_file:
        for chrom in chromsomes_in_vcf.stdout.split('\n'):
            chrom_size = chrom_sizes.get(chrom, default_chromsize)
            if chrom is None:
                raise RuntimeError("Chromsome {} in vcf file but not in "
                                   "chromsome sizes list".format(chrom))
            chrom_region = "{seq}:{start:d}-{end:d}".format(
                seq=chrom, start=1, end=chrom_size
            )
            callset = allel.read_vcf(vcf_file, region=chrom_region)
            if 'calldata/GT' not in callset:
                continue
            archaic_sample_idx = [callset['samples'].tolist().index(i)
                                  for i in archaic_sample_list]
            logging.debug("Archaic sample indicies: {}"
                          .format(archaic_sample_idx))
            modern_sample_idx = [callset['samples'].tolist().index(i)
                                 for i in modern_sample_list]
            logging.debug("Modern sample indicies: {}"
                          .format(modern_sample_idx))
            genotype_array = allel.GenotypeArray(callset['calldata/GT'])
            variant_pos = callset['variants/POS']
            for window_coords in generate_windows(0, chrom_size,
                                                  args.window_size,
                                                  args.step_size):
                region = "{seq}:{start:d}-{end:d}".format(
                    seq=chrom, start=window_coords[0], end=window_coords[1]
                )
                logging.debug("window region: {}".format(region))
                window_length = window_coords[1] - window_coords[0]
                variants_in_region = numpy.where(numpy.logical_and(
                    variant_pos > window_coords[0],
                    variant_pos <= window_coords[1]))
                window_genotype_array = genotype_array.subset(
                    sel0=variants_in_region[0])
                allele_counts = (window_genotype_array
                                 .count_alleles_subpops(subpops={
                                     'archaic': archaic_sample_idx,
                                     'modern': modern_sample_idx}))
                # Find variants with at least one non-reference allele call
                archaic_variant_sites = allele_counts['archaic'].is_variant()
                # Find segregating variants (where more than one allele is
                # observed)
                modern_segregating_sites = (allele_counts['modern']
                                            .is_segregating())

                # Account for masked bases in window
                informative_sites = (archaic_variant_sites &
                                     modern_segregating_sites)
                informative_site_frequency = (sum(informative_sites) /
                                              window_length)
                window = Window(
                    seqid=chrom,
                    start=window_coords[0],
                    end=window_coords[1],
                    informative_site_frequency=informative_site_frequency)
                window_isf_int = int(round(
                    window.informative_site_frequency
                    * args.window_size))
                threshold_int = int(round(args.frequency_threshold *
                                          args.window_size))
                logging.debug("Window: {}".format(region))
                logging.debug("Number of archic variant sites: {}".format(
                              sum(archaic_variant_sites)))
                logging.debug("Number of modern segregating sites: {}".format(
                              sum(modern_segregating_sites)))
                logging.debug("Number of informative sites: {}"
                              .format(sum(informative_sites)))
                # windows_file.write("{}\n".format(window.to_bed()))

                # For each modern haplotype, compare to each archaic haplotype
                modern_haplotypes = allel.HaplotypeArray(
                    window_genotype_array[:, modern_sample_idx]
                    .flatten().reshape(-1, len(modern_sample_idx) * 2))
                archic_haplotypes = allel.HaplotypeArray(
                    window_genotype_array[:, archaic_sample_idx]
                    .flatten().reshape(-1, len(archaic_sample_idx) * 2))
                for hap_idx, modern_haplotype in \
                        enumerate(modern_haplotypes.T):
                    modern_sample = modern_sample_list[
                        int(numpy.trunc(hap_idx / 2))]
                    modern_haplotype_id = "{sample}:{idx}".format(
                        sample=modern_sample,
                        idx=hap_idx % 2)
                    max_match_pct = 0
                    for archic_haplotype in archic_haplotypes.T:
                        num_informative_sites = numpy.sum(informative_sites)
                        if num_informative_sites == 0:
                            match_pct = 0
                        else:
                            match_pct = (
                                numpy.sum(informative_sites &
                                          (modern_haplotype ==
                                           archic_haplotype)) /
                                numpy.sum(informative_sites))
                        if match_pct > max_match_pct:
                            max_match_pct = match_pct
                    if args.match_pct_database:
                        pvalue = max_match_pct_pvalue(
                            window_isf_int,
                            population=(sample_populations[modern_sample]
                                        .population),
                            max_match_pct=max_match_pct,
                            dbconn=dbconn,
                            threshold=threshold_int)
                        print("{chr}\t{start}\t{end}\t{isf:.6f}\t{isfd:d}\t"
                              "{haplotype}\t{population}\t{match_pct}\t"
                              "{pvalue}".format(
                                  chr=window.seqid,
                                  start=window.start,
                                  end=window.end,
                                  isf=window.informative_site_frequency,
                                  isfd=window_isf_int,
                                  haplotype=modern_haplotype_id,
                                  population=(
                                      sample_populations[modern_sample]
                                      .population),
                                  match_pct=max_match_pct,
                                  pvalue=pvalue))
                    else:
                        print("{isfd:d}\t{population}\t{match_pct}".format(
                            isfd=window_isf_int,
                            population=sample_populations[modern_sample]
                            .population,
                            match_pct=max_match_pct))
                sys.stdout.flush()
                logging.debug(
                    "windows_within_isf_threshold.cache_info(): {}"
                    .format(windows_within_isf_threshold.cache_info()))
                logging.debug(
                    "max_match_pct_pvalue.cache_info(): {}"
                    .format(max_match_pct_pvalue.cache_info()))
        sys.stdout.flush()


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
