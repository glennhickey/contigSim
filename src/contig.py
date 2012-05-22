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
contigs have an extra edge because of implicit telomere bases.  the operations
in the member functions don't change the object (return new ones instead).
New telomeres aren't created when cutting a linear contig (done outside)

"""
class Contig(object):
    def __init__(self, size):
        self.size = size

    def setDead(self, dead=True):
        self.dead = dead

    def isDead(self):
        if hasattr(self, 'dead'):
            return self.dead
        return False

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

    # cut (remove edge) at position, returning a left contig and a right contig
    # (original not touched)
    def cut(self, position):
        assert position < self.size
        l = LinearContig(position)
        r = LinearContig(self.size - position - 1)
        return (l,r)

    # get the contig in reverse orientation
    def reverse(self):
        return copy.deepcopy(self)

    # add edge between two endpoints
    # return the new circular contig (current contig is unchange)
    def circularize(self):
        return CircularContig(self.size + 1)

    # stick another linear contig to the left (new edge)
    # forward is the direction of the other contig
    def joinToLeft(self, other, forward=True):
        assert type(other) == LinearContig
        return LinearContig(self.size + other.size + 1)

    # stick another linear contig to the right (new edge)
    # forward is the direction of the other contig
    def joinToRight(self, other, forward=True):
        assert type(other) == LinearContig
        return LinearContig(self.size + other.size + 1)



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

    # chop a contig into two circles
    def cut(self, pos1, pos2):
        left = CircularContig(self.size - math.fabs(pos1 - pos2))
        right = CircularContig(math.fabs(pos1 - pos2))
        return (left, right)
    
    # return a linearized version of the circular contig by cutting at position
    # (removing edge)
    def linearize(self, position = 0):
        return LinearContig(self.size - 1)

    # join two circles with a dcj operation
    def join(self, other, pos1 = 0, pos2 = 0, forward = True):
        return CircularContig(self.size + other.size)
    

    
    
        

    
