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
from .classmodule import Window
from .funcmodule import get_samplename_list
from .funcmodule import get_chrom_sizes
from .funcmodule import generate_windows

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
    input_parser = argparse.ArgumentParser(add_help=False)
    input_group = input_parser.add_argument_group('Variant input files')
    input_group.add_argument("--vcf-file",
                             required=True,
                             help="VCF file to parse",
                             metavar="FILE")
    input_group.add_argument("--archaic-samples",
                             required=True,
                             help="File listing archaic samples, one per line",
                             metavar="FILE")
    input_group.add_argument("--modern-samples",
                             required=True,
                             help="File listing modern samples to analyze, "
                             "one per line",
                             metavar="FILE")
    input_group.add_argument("--chrom-sizes",
                             required=True,
                             help="Tab delimited file with seqid\tlength",
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
                              default=0.0005,
                              type=float,
                              metavar="THRESHOLD")

    # Window size options
    window_parser = argparse.ArgumentParser(add_help=False)
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
    parent_parser.add_argument("-v", "--verbose",
                               help="increase output verbosity",
                               action="store_true")
    parent_parser.add_argument("--version", action="version",
                               version="%(prog)s " + __version__)

    subparsers = parent_parser.add_subparsers(
        help="commands")
    parser_max_match_pct = subparsers.add_parser(
        'max-match-pct',
        parents=[input_parser, window_parser])
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


def max_match_pct_pvalue(informative_site_frequency, haplotype,
                         max_match_pct, dbconn, threshold):
    '''Calculate emperical pvalue from windows in database'''
    lower_threshold = informative_site_frequency - threshold
    upper_threshold = informative_site_frequency + threshold
    c = dbconn.cursor()
    c.execute('''SELECT max_match_pct from MAX_MATCH_PCT
              where informative_site_frequency >= ?
              and informative_site_frequency <= ?''',
              (lower_threshold, upper_threshold))
    null_distribution = numpy.array(c.fetchall())
    pvalue = float(sum(null_distribution >= max_match_pct)
                   / len(null_distribution))
    return pvalue


def max_match_pct(args):

    # List of archic sample names
    archaic_sample_list = get_samplename_list(args.archaic_samples)
    logging.debug("Archaic sample list: {}".format(archaic_sample_list))

    modern_sample_list = get_samplename_list(args.modern_samples)
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

    chromsomes_in_vcf = subprocess.run(['tabix', '-l', args.vcf_file],
                                       stdout=subprocess.PIPE,
                                       universal_newlines=True)
    # logging.debug(chromsomes_in_vcf.stdout)
#    with open(args.windows_file, 'w') as windows_file:
    for chrom in chromsomes_in_vcf.stdout.split('\n'):
        chrom_size = chrom_sizes.get(chrom, default_chromsize)
        if chrom is None:
            raise RuntimeError("Chromsome {} in vcf file but not in chromsome "
                               "sizes list".format(chrom))
        chrom_region = "{seq}:{start:d}-{end:d}".format(
            seq=chrom, start=1, end=chrom_size
        )
        callset = allel.read_vcf(args.vcf_file, region=chrom_region)
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
            # TODO This should be where one or more archaic haplotypes has derived
            archaic_variant_sites = allele_counts['archaic'].is_variant()
            # TODO This should be where modern pop varies
            modern_variant_sites = allele_counts['modern'].is_variant()

            # TODO Account for masked bases in window
            informative_sites = (archaic_variant_sites &
                                 modern_variant_sites)
            informative_site_frequency = (sum(informative_sites) /
                                          window_length)
            window = Window(
                seqid=chrom,
                start=window_coords[0],
                end=window_coords[1],
                informative_site_frequency=informative_site_frequency)
            logging.debug("Window: {}".format(region))
            logging.debug("Number of archic variant sites: {}".format(
                          sum(archaic_variant_sites)))
            logging.debug("Number of moden variant sites: {}".format(
                          sum(modern_variant_sites)))
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
                    # TODO Limit to informative_sites
                    match_pct = (numpy.sum(
                        modern_haplotype == archic_haplotype) /
                        window_length)
                    if match_pct > max_match_pct:
                        max_match_pct = match_pct
                if args.match_pct_database:
                    pvalue = max_match_pct_pvalue(
                        window.informative_site_frequency,
                        haplotype=modern_haplotype_id,
                        max_match_pct=max_match_pct,
                        dbconn=dbconn,
                        threshold=args.frequency_threshold)
                    print("{window:.6f}\t{haplotype}\t{match_pct}\t{pvalue}"
                          .format(window=window.informative_site_frequency,
                                  haplotype=modern_haplotype_id,
                                  match_pct=max_match_pct,
                                  pvalue=pvalue))
                else:
                    print("{window:.6f}\t{haplotype}\t{match_pct}".format(
                        window=window.informative_site_frequency,
                        haplotype=modern_haplotype_id,
                        match_pct=max_match_pct))
    sys.stdout.flush()


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
