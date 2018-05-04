#!/usr/bin/env python
# -*- coding: utf-8 -*-


class Window():
    """docstring for [object Object]."""

    def __init__(self, seqid, start, end, informative_sites_count=None):
        if start < 0:
            raise AttributeError("Window start must be greater than 0: {} < 0"
                                 .format(start))
        if start >= end:
            raise AttributeError("Window start must be less than window end: "
                                 "{} >= {}".format(start, end))
        self.seqid = seqid
        self.start = start
        self.end = end
        self._region_string = "{seqid}:{start}-{end}".format(
            seqid=self.seqid, start=self.start, end=self.end)
        self._size = self.end - self.start
        self._informative_sites_count = informative_sites_count
        self._informative_sites_frequency = None

    @property
    def size(self):
        return self._size

    @property
    def informative_sites_frequency(self):
        if (self._informative_sites_frequency is None
                and self._informative_sites_count is not None):
            self._informative_sites_frequency = (
                self._informative_sites_count / self.size
            )
        return self._informative_sites_frequency

    @property
    def informative_sites_count(self):
        return self._informative_sites_count

    @informative_sites_count.setter
    def informative_sites_count(self, value):
        if value > self.size:
            raise AttributeError("informative_sites_count cannot be greater "
                                 "that the size of the window: {} > {}".format(
                                     value, self.size
                                 ))
        self._informative_sites_count = value
        self._informative_sites_frequency = (
            self._informative_sites_count / self.size
        )

    @property
    def region_string(self):
        return self._region_string

    def to_bed(self):
        return("{seqid}\t{start}\t{end}\t{isc}"
               .format(seqid=self.seqid,
                       start=self.start,
                       end=self.end,
                       isc=self.isc))
