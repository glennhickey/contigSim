#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import copy
import math

""" a contig is an interval of DNA that can either be circular or linear.  a
contig can be decomposed into an alternating walk of bases and adjacency
edges. in the current implementation we only keep track of the number of
adjacency edges in the contig, and whether its circular or linear. linear
contigs have an extra edge because of implicit telomere bases.  

"""
class Contig(object):
    def __init__(self, size):
        self.size = size

class LinearContig(Contig):
    def __init__(self, size = 1):
         super(LinearContig, self).__init__(size)

    def isLinear(self):
        return True

    def isCircular(self):
        return False

    # contig has form o-x-x-x-o (where -'s are adjacencies, x's are bases
    # and o's are telomeres).  size reflects the number of -'s.  so the
    # number of bases is size - 1
    def numBases(self):
        if self.size <= 1:
            return 0
        else:
            return self.size - 1

    # break the contig at its positionth adjacency.  return a new contig
    # of everything to the left (and keep everything to the right).  two
    # new telomeres are implicity created
    def cutLeftOff(self, position):
        assert position < self.size
        self.size -= position
        return LinearContig(position + 1)

    # break the contig at its positionth adjacency.  return a new contig
    # of everything to the right (and keep everything to the left).  two
    # new telomeres are implicity created
    def cutRightOff(self, position):
        assert position < self.size
        newSize = self.size - position
        self.size = position + 1
        return LinearContig(newSize)

    # pop off the telomeres and close the linear contig into a circle
    # return the new circular contig (current contig is unchange)
    def circularize(self):
        return CircularContig(self.size - 1)

    # stick another linear contig to the left
    def joinToLeft(self, other):
        assert type(other) == LinearContig
        self.size += other.size - 1

    # stick another linear contig to the right
    # same as above since without actual bases its symmetric
    def joinToRight(self, other):
        assert type(other) == LinearContig
        self.size += other.size - 1


class CircularContig(Contig):
    def __init__(self, size = 0):
         super(CircularContig, self).__init__(size)

    def isLinear(self):
        return False

    def isCircular(self):
        return True

    # contig has form o-x-x-x-o (where -'s are adjacencies, x's are bases
    # and o's are telomeres).  size reflects the number of -'s.  so the
    # number of bases is size - 1
    def numBases(self):
        if self.size <= 0:
            return 0
        else:
            return self.size

    # chop off a segment and make a new circular contig out of it
    def cutOffCircle(self, pos1, pos2):
        self.size -= math.fabs(pos1 - pos2)
        return CircularContig(math.fabs(pos1 - pos2))
    
    # return a linearized version of the circular contig by cutting at position
    def linearize(self, position = 0):
        return LinearContig(self.size + 1)

    # cut cricle at pos1 (AB) and other circle at pos2 (CD)
    # join them to form new adjacencies AC and BD
    def joinOval(self, other, pos1 = 0, pos2 = 0):
        assert type(other) == CircularContig
        self.size += other.size
        
    # cut cricle at pos1 (AD) and other circle at pos2 (BC)
    # join them to form new adjacencies AC and BD
    def joinEight(self, other, pos1 = 0, pos2 = 0):
        assert type(other) == CircularContig
        self.size += other.size


    
    
        

    
