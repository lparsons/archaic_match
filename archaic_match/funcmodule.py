#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from collections import namedtuple
from collections import defaultdict
from .classmodule import Window


def get_samplename_list(query_populations, sample_populations):
    '''Returns sample names for given population'''
    sample_list = list()
    for sample in sample_populations.values():
        if sample.population in query_populations:
            sample_list.append(sample.name)
    return sample_list


def get_chrom_sizes(filename):
    '''Return dictionary with seqid as key and length as value'''
    d = defaultdict()
    with open(filename) as fin:
        rows = (line.split('\t') for line in fin)
        d = {row[0]: int(row[1]) for row in rows}
    return d


def tsv_to_dictionary(filename):
    '''Return dictionary with first column as key and second as value'''
    d = defaultdict()
    with open(filename) as fin:
        rows = (line.rstrip('\n').split('\t') for line in fin)
        d = {row[0]: row[1] for row in rows}
    return d


Sample = namedtuple('Sample', ['name', 'population', 'superpopulation'])


def get_sample_populations(filename):
    '''Return dictionary {sample: (name, population, superpopulation)}'''
    with open(filename) as fin:
        rows = (line.split('\t') for line in fin)
        d = {row[0]: Sample(name=row[0], population=row[1],
                            superpopulation=row[2].strip()) for row in rows}
    return d


def generate_windows(seqid, start, end, size, step):
    '''Generate windows stepping along chromosome'''
    for winstart in range(start, end, step):
        winend = min(winstart + size, end)
        yield Window(seqid, winstart, winend)


def get_informative_sites(allele_counts, method):
    """
    Determine which sites are considered informative
    subpopulation allele counts for "archaic" and "modern"

    Parameters
    ----------
    allele_counts : dict (string -> AlleleCountsArray)
    method : str
        One of `derived_in_archaic` or `derived_in_archaic_or_modern`

    Returns
    -------
    out : ndarray, bool, shape (n_variants,)
    """

    if method == "derived_in_archaic":
        informative_sites = allele_counts['archaic'].is_variant()
        logging.debug("Informative site method: {}".format(method))
        logging.debug("Number of archic variant sites: {}".format(
            sum(informative_sites)))
        logging.debug("Number of informative sites: {}".format(
            sum(informative_sites)))
    elif method == "derived_in_archaic_or_modern":
        archaic_variant_sites = allele_counts['archaic'].is_variant()
        modern_variant_sites = allele_counts['modern'].is_variant()
        informative_sites = (archaic_variant_sites | modern_variant_sites)
        logging.debug("Informative site method: {}".format(method))
        logging.debug("Number of archic variant sites: {}".format(
            sum(archaic_variant_sites)))
        logging.debug("Number of modern variant sites: {}".format(
            sum(modern_variant_sites)))
        logging.debug("Number of informative sites: {}".format(
            sum(informative_sites)))
    else:
        raise RuntimeError("Informative site method '{}' is not supported",
                           method)

    return informative_sites
