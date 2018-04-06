#!/usr/bin/env python
# -*- coding: utf-8 -*-


def my_function(text_to_display):
    print('text from my_function :: {}'.format(text_to_display))


def get_samplename_list(filename):
    '''Return lines in file as a list'''
    with open(filename) as f:
        sample_list = f.readlines()
        sample_list = [x.strip() for x in sample_list]
    return sample_list


def get_chrom_sizes(filename):
    '''Return dictionary with seqid as key and length as value'''
    with open(filename) as fin:
        rows = (line.split('\t') for line in fin)
        d = {row[0]: int(row[1]) for row in rows}
    return d


def generate_windows(start, end, size, step):
    '''Return list of windows as tuples'''
    for winstart in range(start, end, step):
        winend = min(winstart + size, end)
        yield (winstart, winend)
