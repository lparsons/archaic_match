#!/usr/bin/env python
# -*- coding: utf-8 -*-


class MyClass():
    def __init__(self, name):
        self.name = name

    def say_name(self):
        print('name is {}'.format(self.name))


class Window():
    """docstring for [object Object]."""

    def __init__(self, seqid, start, end, informative_site_frequency):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.informative_site_frequency = informative_site_frequency

    def to_bed(self):
        return("{seqid}\t{start}\t{end}\t{isf}"
               .format(seqid=self.seqid,
                       start=self.start,
                       end=self.end,
                       isf=self.informative_site_frequency))
