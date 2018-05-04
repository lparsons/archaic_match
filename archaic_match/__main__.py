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
import pandas
from collections import namedtuple
from collections import defaultdict
from .classmodule import Window
from .funcmodule import get_samplename_list
from .funcmodule import get_chrom_sizes
from .funcmodule import generate_windows
from .funcmodule import get_sample_populations
from .funcmodule import get_informative_sites


__version__ = "0.2"
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
    pvalue_group.add_argument(
        "--overlap-regions",
        required=False,
        help="Report number of basepairs and number of informative sites "
        "of each window that overlap these regions",
        metavar="BED_FILE"
    )

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

    # Algoritmic options
    algorithm_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    algorithm_group = algorithm_parser.add_argument_group('Algorithm options')
    algorithm_group.add_argument("--informative-site-method",
                                 help=("Specify method used to define "
                                       "sites"),
                                 default="derived_in_archaic",
                                 choices=["derived_in_archaic",
                                          "dervived_in_archaic_or_modern"])

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
        parents=[input_parser, window_parser, common_options],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Calculate match percent values for modern haplotypes")
    parser_max_match_pct.set_defaults(func=max_match_pct)

    # Build database subcommand
    build_db_parser = subparsers.add_parser(
        'build-db',
        parents=[common_options],
        help="Build database from match percent count files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    build_db_parser.set_defaults(func=build_db)

    build_req_args = build_db_parser.add_argument_group("required arguments")
    build_req_args.add_argument(
        "--match-pct-count",
        help="Match percent count file(s) to load",
        required=True,
        nargs='+',
        metavar='FILE'
    )
    build_req_args.add_argument(
        "--db",
        help="Database file",
        required=True,
        metavar='DB_FILE'
    )

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


def build_db(args):
    logging.debug("Creating database in file: '{}'".format(args.db))

    import csv
    count_filelist = list()
    for count_file_pattern in args.match_pct_count:
        count_filelist.extend(glob.glob(count_file_pattern))
    logging.debug("Input file list:\n{}".format(count_filelist))

    logging.debug("Creating match_pct_counts table")
    import sqlite3
    dbconn = sqlite3.connect(args.db)
    c = dbconn.cursor()
    c.execute('DROP TABLE IF EXISTS match_pct_counts')
    c.execute('CREATE TABLE match_pct_counts('
              'informative_site_frequency INTEGER, '
              'population TEXT, '
              'max_match_pct REAL, '
              'count INTEGER)')

    logging.debug("Creating tmp table")
    c.execute('DROP TABLE IF EXISTS tmp')
    c.execute('CREATE TABLE tmp('
              'informative_site_frequency INTEGER, '
              'population TEXT, '
              'max_match_pct REAL, '
              'count INTEGER)')
    for count_file in count_filelist:
        logging.debug("Loading file '{}'".format(count_file))
        tsvData = csv.reader(open(count_file, "rU"), delimiter="\t")
        divData = chunks(tsvData)  # divide into 10000 rows each
        for chunk in divData:
            c.executemany('INSERT OR IGNORE INTO tmp '
                          '(informative_site_frequency, population, '
                          'max_match_pct, count) '
                          'VALUES (?,?,?,?)', chunk)
    dbconn.commit()

    logging.debug("Summarize temp data into match_pct_counts from tmp table")
    c.execute('insert into match_pct_counts '
              '(informative_site_frequency, population, '
              'max_match_pct, count) '
              'select distinct informative_site_frequency, '
              'population, max_match_pct, sum(count) from tmp '
              'group by informative_site_frequency, population, '
              'max_match_pct')
    dbconn.commit()
    c.execute('DROP TABLE tmp')
    logging.debug("Creating indexes")
    c.execute('CREATE INDEX IF NOT EXISTS isf_idx ON match_pct_counts '
              '(informative_site_frequency)')
    c.execute('CREATE INDEX IF NOT EXISTS mmp_idx ON match_pct_counts '
              '(max_match_pct)')
    c.execute('CREATE INDEX IF NOT EXISTS pop_idx ON match_pct_counts '
              '(population)')
    c.execute('CREATE INDEX IF NOT EXISTS isf_pop_idx ON match_pct_counts '
              '(informative_site_frequency, population)')
    c.execute('CREATE INDEX IF NOT EXISTS isf_pop_mmpct_idx ON '
              'match_pct_counts (informative_site_frequency, population, '
              'max_match_pct)')
    logging.debug("Analyzing")
    c.execute('ANALYZE')
    logging.debug("Vacuuming")
    c.execute('VACUUM')
    dbconn.commit()


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
    c.execute('''SELECT sum(count) from match_pct_counts
              where informative_site_frequency >= ?
              and informative_site_frequency <= ?
              and population = ?
              and max_match_pct >= ?''',
              (lower_threshold, upper_threshold, population, max_match_pct))
    s = c.fetchone()[0]
    if s is None:
        s = 0
    s = float(s)
    pvalue = (s + 1) / (n + 1)
    return (pvalue, int(n))


@functools.lru_cache(maxsize=None)
def windows_within_isf_threshold(informative_site_frequency, population,
                                 threshold, dbconn):
    '''Return number of windows for a given population within the informative
    site frequency threshold'''
    lower_threshold = informative_site_frequency - threshold
    upper_threshold = informative_site_frequency + threshold
    c = dbconn.cursor()
    c.execute('''SELECT sum(count) from match_pct_counts
              where informative_site_frequency >= ?
              and informative_site_frequency <= ?
              and population = ?''',
              (lower_threshold, upper_threshold, population))
    n = c.fetchone()[0]
    if n is None:
        n = 0
    n = float(n)
    return n


match_pct_window_pval = namedtuple(
    'match_pct_window_pval',
    ['chr', 'start', 'end', 'isf', 'isfd', 'haplotype', 'population',
     'match_pct', 'pvalue', 'matching_windows', 'overlap_bp',
     'overlap_informative_sites'])


match_pct_window = namedtuple(
    'match_pct_window', ['isfd', 'population', 'match_pct'])


# match_pct_window_counts = namedtuple(
#     'match_pct_window_counts',
#     ['isfd', 'population', 'match_pct', 'count'])


def calc_match_pct(informative_sites, archaic_haplotypes, modern_haplotype):
    '''Cacluate maximum match percent of modern haplotype to archaic
    haplotypes'''
    archaic_is_derived = archaic_haplotypes.count_alleles().is_variant()
    modern_is_derived = modern_haplotype > 0
    match_pct = numpy.sum(
        (informative_sites & archaic_is_derived & modern_is_derived)
        / numpy.sum(informative_sites))
    return match_pct


def calc_window_haplotype_match_pcts(
        vcf_file, chrom_sizes,
        archaic_sample_list, modern_sample_list,
        window_size, step_size, frequency_threshold,
        informative_site_method,
        dbconn=None, sample_populations=None, overlap_regions=None):
    '''Generate match pct for each window for each modern haplotype'''
    # TODO Use pysam here
    # http://pysam.readthedocs.io/en/latest/api.html#pysam.TabixFile
    chromsomes_in_vcf = (subprocess.run(['tabix', '-l', vcf_file],
                                        stdout=subprocess.PIPE,
                                        universal_newlines=True)
                         .stdout.split('\n'))
    for chrom in chromsomes_in_vcf:
        chrom_size = chrom_sizes[chrom]
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

        logging.debug(overlap_regions)
        if (overlap_regions is not None) and (not overlap_regions.empty):
            chrom_overlap_regions_bool = overlap_regions['chr'] == str(chrom)
            chrom_overlap_regions = overlap_regions[chrom_overlap_regions_bool]
        else:
            chrom_overlap_regions = overlap_regions
        for window_coords in generate_windows(0, chrom_size,
                                              window_size,
                                              step_size):
            region = "{seq}:{start:d}-{end:d}".format(
                seq=chrom, start=window_coords[0], end=window_coords[1]
            )
            logging.debug("window region: {}".format(region))
            window_length = window_coords[1] - window_coords[0]
            variants_in_region = numpy.where(numpy.logical_and(
                variant_pos > window_coords[0],
                variant_pos <= window_coords[1]))
            variant_pos_in_region = variant_pos[variants_in_region]
            window_genotype_array = genotype_array.subset(
                sel0=variants_in_region[0])
            allele_counts = (window_genotype_array
                             .count_alleles_subpops(subpops={
                                 'archaic': archaic_sample_idx,
                                 'modern': modern_sample_idx}))
            informative_sites = get_informative_sites(
                allele_counts, informative_site_method)
            informative_site_frequency = (sum(informative_sites)
                                          / window_length)
            window = Window(
                seqid=chrom,
                start=window_coords[0],
                end=window_coords[1],
                informative_site_frequency=informative_site_frequency)
            window_isf_int = int(round(
                window.informative_site_frequency * window_size))
            threshold_int = int(round(frequency_threshold
                                      * window_size))

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

                max_match_pct = calc_match_pct(
                    informative_sites, archic_haplotypes, modern_haplotype)

                if dbconn:
                    (pvalue, matching_windows) = max_match_pct_pvalue(
                        window_isf_int,
                        population=(sample_populations[modern_sample]
                                    .population),
                        max_match_pct=max_match_pct,
                        dbconn=dbconn,
                        threshold=threshold_int)
                    # TODO Calculate overlap with optional region file
                    (overlap_bp, overlap_informative_sites) = \
                        calculate_overlap(
                            chrom_overlap_regions, window, modern_haplotype_id,
                            variant_pos_in_region[informative_sites])
                    yield match_pct_window_pval(
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
                        pvalue=pvalue,
                        matching_windows=matching_windows,
                        overlap_bp=overlap_bp,
                        overlap_informative_sites=overlap_informative_sites)
                else:
                    yield match_pct_window(
                        isfd=window_isf_int,
                        population=sample_populations[modern_sample]
                        .population,
                        match_pct=max_match_pct)
        logging.debug(
            "windows_within_isf_threshold.cache_info(): {}"
            .format(windows_within_isf_threshold.cache_info()))
        logging.debug(
            "max_match_pct_pvalue.cache_info(): {}"
            .format(max_match_pct_pvalue.cache_info()))


def calculate_overlap(chrom_overlap_regions, window, modern_haplotype_id,
                      informative_site_positions):
    overlapping_bp = 0
    overlapping_informative_sites = list()
    if not chrom_overlap_regions.empty:
        sample_chrom_overlap_regions = (
            chrom_overlap_regions
            [chrom_overlap_regions['sample'] == modern_haplotype_id])
        if not sample_chrom_overlap_regions.empty:
            logging.debug("Overlap regions in chrom:\n{}".format(
                chrom_overlap_regions))
            logging.debug("Window: {}".format(window.start))
            overlapping_regions = sample_chrom_overlap_regions[
                (chrom_overlap_regions['start'] <= window.end)
                & (chrom_overlap_regions['end'] >= window.start)]
            logging.debug("Overlapping regions for window:\n{}".format(
                overlapping_regions))
            overlapping_bp = 0
            for index, region in overlapping_regions.iterrows():
                logging.debug(region)
                overlap = (min(region['end'], window.end)
                           - max(region['start'], window.start))
                overlapping_bp += overlap
            logging.debug("Informative site positions: {}".format(
                informative_site_positions))
            informative_site_index = allel.SortedIndex(
                informative_site_positions)
            overlapping_informative_sites = (
                informative_site_index.intersect_ranges(
                    starts=overlapping_regions['start'],
                    stops=overlapping_regions['end']))
    return(overlapping_bp, len(overlapping_informative_sites))


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
        chrom_sizes = defaultdict(lambda: int(args.chrom_sizes))
    else:
        chrom_sizes = get_chrom_sizes(args.chrom_sizes)
    logging.debug("Chromosome sizes: {}".format(chrom_sizes))

    dbconn = None
    if args.match_pct_database:
        import sqlite3
        dbconn = sqlite3.connect(args.match_pct_database)

    vcf_filelist = list()
    for vcf_file_pattern in args.vcf:
        vcf_filelist.extend(glob.glob(vcf_file_pattern))
    logging.debug("VCF File list:\n{}".format(vcf_filelist))
    match_pct_window_counts = defaultdict(int)
    overlap_regions = None
    if dbconn:
        if args.overlap_regions:
            fields = match_pct_window_pval._fields
        else:
            fields = match_pct_window_pval._fields[:-2]
        print("\t".join(fields))
        region_colnames = ['chr', 'start', 'end', 'sample']
        region_dtypes = {'chr': object,
                         'start': int,
                         'end': int,
                         'sample': object}
        if args.overlap_regions:
            overlap_regions = pandas.read_table(args.overlap_regions,
                                                header=None,
                                                names=region_colnames,
                                                dtype=region_dtypes)
        else:
            overlap_regions = pandas.DataFrame(columns=region_colnames)
    for vcf_file in vcf_filelist:
        window_haplotype_match_pcts = calc_window_haplotype_match_pcts(
            vcf_file=vcf_file, chrom_sizes=chrom_sizes,
            archaic_sample_list=archaic_sample_list,
            modern_sample_list=modern_sample_list,
            window_size=args.window_size,
            step_size=args.step_size,
            frequency_threshold=args.frequency_threshold,
            dbconn=dbconn,
            sample_populations=sample_populations,
            overlap_regions=overlap_regions)
        if dbconn:
            for window_haplotype_match_pct in window_haplotype_match_pcts:
                sys.stdout.write(
                    "{chr}\t{start}\t{end}\t{isf:.6f}\t{isfd:d}\t"
                    "{haplotype}\t{population}\t{match_pct}\t"
                    "{pvalue}\t{matching_windows}"
                    .format(**window_haplotype_match_pct._asdict()))
                if args.overlap_regions:
                    sys.stdout.write(
                        "\t{overlap_bp}\t{overlap_informative_sites}"
                        .format(**window_haplotype_match_pct._asdict()))
                sys.stdout.write("\n")
        else:
            for window_haplotype_match_pct in window_haplotype_match_pcts:
                match_pct_window_counts[window_haplotype_match_pct] += 1
    if not dbconn:
        for w in match_pct_window_counts.keys():
            print("{isfd}\t{population}\t{match_pct}\t{count}"
                  .format(**w._asdict(), count=match_pct_window_counts[w]))
    sys.stdout.flush()


def chunks(reader, chunksize=10000):
    """
    Chunk generator. Take a CSV `reader` and yield
    `chunksize` sized slices.
    """
    chunk = []
    for i, line in enumerate(reader):
        if (i % chunksize == 0 and i > 0):
            yield chunk
            del chunk[:]
        chunk.append(line)
    yield chunk


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
