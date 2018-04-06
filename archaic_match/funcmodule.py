#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import namedtuple


def my_function(text_to_display):
    print('text from my_function :: {}'.format(text_to_display))


def get_samplename_list(query_populations, sample_populations):
    '''Returns sample names for given population'''
    sample_list = list()
    for sample in sample_populations.values():
        if sample.population in query_populations:
            sample_list.append(sample.name)
    return sample_list


def get_chrom_sizes(filename):
    '''Return dictionary with seqid as key and length as value'''
    with open(filename) as fin:
        rows = (line.split('\t') for line in fin)
        d = {row[0]: int(row[1]) for row in rows}
    return d


Sample = namedtuple('Sample', ['name', 'population', 'superpopulation'])


def get_sample_populations(filename):
    '''Return dictionary {sample: (name, population, superpopulation)}'''
    with open(filename) as fin:
        rows = (line.split('\t') for line in fin)
        d = {row[0]: Sample(name=row[0], population=row[1],
                            superpopulation=row[2].strip()) for row in rows}
    return d


def generate_windows(start, end, size, step):
    '''Return list of windows as tuples'''
    for winstart in range(start, end, step):
        winend = min(winstart + size, end)
        yield (winstart, winend)
