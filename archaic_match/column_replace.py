#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" column_replace

Replace column values using lookup files

This is useful to rename chromsomes in a bed file, for example.
"""

import argparse
import logging
import glob
import pandas
import sys
from .funcmodule import tsv_to_dictionary

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2018, Lance Parsons"
__license__ = "MIT https://opensource.org/licenses/MIT"


def main():

    # Top level parser
    parser = argparse.ArgumentParser(
        description="Replace column values using a lookup file",
        epilog="As an alternative to the commandline, params can be placed in "
        "a file, one per line, and specified on the commandline like "
        "'%(prog)s @params.conf'.",
        fromfile_prefix_chars='@',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True)

    parser.add_argument(
        '-d', '--dictionary',
        required=True,
        help="Dictionary tsv file, with column for current value and "
        "replacement value (may also be a quoted glob expression)",
        metavar="DICTIONARY_FILE"
    )

    parser.add_argument(
        '-c', '--column',
        help="Column where values will be replaced",
        default=1,
        type=int,
        metavar="COLUMN_NUMBER"
    )

    parser.add_argument(
        'input_files',
        nargs="+",
        help="TSV files in which to replace values",
        metavar="TSV_FILE")

    parser.add_argument(
        "-v", "--verbose",
        help="increase output verbosity",
        action="store_true")
    parser.add_argument(
        "--version", action="version",
        version="%(prog)s " + __version__)

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    input_file_list = list()
    for input_file_pattern in args.input_files:
        input_file_list.extend(glob.glob(input_file_pattern))

    logging.debug("Replacing values in column {} of the following files:\n{}"
                  .format(args.column, input_file_list))
    logging.debug("Dictionary file: '{}'".format(args.dictionary))
    logging.debug("Values matching something column 1 of the dictionary will "
                  "be replaced with the corrsponding value in column 2 of the "
                  "dictionary")

    replace_dictionary = tsv_to_dictionary(args.dictionary)
    logging.debug(replace_dictionary)

    col_num = args.column - 1
    for input_file in input_file_list:
        try:
            reader = pandas.read_table(
                input_file, header=None, chunksize=1000,
                converters={col_num: lambda x: replace_dictionary.get(x, x)})
            for df in reader:
                df.to_csv(sys.stdout, sep="\t", header=False, index=False)
        except pandas.io.common.EmptyDataError:
            logging.debug("'{}' is empty and was skipped.".format(input_file))
