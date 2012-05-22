#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import copy
import random
import math
from heapq import heappush, heappop

from contig import CircularContig
from contig import LinearContig
from contig import Contig
from dcj import dcj
from eventQueue import EventQueue
from sampleTree import SampleTree


# simple rearrangment model with a pool of contigs (circular and linear)
# and a circular garbage contig
class Model(object):
    def __init__(self):
        self.pool = SampleTree()
        self.eventQueue = EventQueue()

    # there are five kinds of rates:
    # N: (fixed) number of bases in the model
    # rll: rate for dcj on the bases in the contig pool
    # rld: rate for dcj where one break is in the pool
    #      and the other rate is in the garbage
    # rdd: both in garbage
    # fl: telomere loss modifier
    # fg: telomere gain modifier
    def setParameters(self, N, rll, rld = 0, rdd = 0, fl = 0, fg = 0):        
        self.eventQueue.reset()
        self.N = N
        self.fl = fl
        self.fg = fg

        if rll > 0:
            self.eventQueue.addEventType(N * rll, self.__llEvent)
        if rld > 0:
            self.eventQueue.addEventType(N * rld, self.__ldEvent)
        if rdd > 0:
            self.eventQueue.addEventType(N * rdd, self.__ddEvent)

    # intitialize the starting state
    # the the contigs will all have the same sizes (modulo rounding)
    # in order to satisfy the input parameters exactly
    def setStartingState(self, garbageSize, numLinear, numCircular):
        assert self.N > garbageSize + numLinear + numCircular
        self.pool = SampleTree()

        numGarbage = 0
        if garbageSize > 0:
            garbage = CircularContig(garbageSize)
            garbage.setDead()
            self.pool.insert(garbage, garbage.size)
            numGarbage = 1
        
        lrat = float(numLinear) / (numLinear + numCircular)
        crat = float(numCircular) / (numLinear + numCircular)
        linearBases = math.floor((self.N - garbageSize) * lrat)
        circularBases = math.ceil((self.N - garbageSize) * crat)
        assert linearBases + circularBases + garbageSize == self.N

        if numLinear > 0:
            linSize = math.floor(linearBases / numLinear)
            extra = linearBases % numLinear
            added = 0
            for i in range(numLinear):
                size = linSize
                if i < extra:
                    size += 1
                # plus 1 since number of adjacencies is 1 + number of bases
                contig = LinearContig(size + 1)
                self.pool.insert(contig, contig.size)
                added += contig.size
            assert added == linearBases + numLinear
            assert self.pool.size() == numLinear + numGarbage
            assert self.pool.weight() == linearBases + numLinear + garbageSize

        if numCircular > 0:
            circSize = math.floor(circularBases / numCircular)
            extra = circularBases % numCircular
            added = 0
            for i in range(numCircular):
                size = circSize
                if i < extra:
                    size += 1
                contig = CircularContig(size)
                self.pool.insert(contig, contig.size)
                added += contig.size
            assert added == circularBases
            assert self.pool.size() == numLinear + numCircular + numGarbage
            assert self.pool.weight() == circularBases + linearBases + \
            numLinear + garbageSize

    # run the simulation for the specified time
    def simulate(self, time):
        self.eventQueue.begin()
        while True:
            nextEvent = self.eventQueue.next(time)
            if nextEvent is not None:
                nextEvent()
            else:
                break
    
    def __llEvent(self):
        if self.pool.size() == 0 or self.pool.weight() == 1:
            return
        
        # draw (and remove) two random adajcenies and their
        #contigs from the pool (only if they are not dead)
        sampleNode1, offset1 = self.pool.uniformSample()
        sampleNode2, offset2 = self.pool.uniformSample()
        c1 = sampleNode1.data
        c2 = sampleNode2.data
        if c1.isDead() == False and c2.isDead() == False:
            self.pool.remove(sampleNode1)
            if c1 is not c2:
                self.pool.remove(sampleNode2)

            gain = None
            # case 1) gain of telomere
            if sampleNode1 is sampleNode2 and offset1 == offset2:
                # correct "not composite check below"
                if c1.isCircular() or (offset1 != 0 and offset1 != c1.size - 1):
                    forward = self.fg > random.random()
                    gain = forward

            # case 2) loss of telomer
            elif c1.isLinear() and c2.isLinear() and \
                    (offset1 == 0 or offset1 == c1.size - 1) and \
                    (offset2 == 0 or offset2 == c2.size - 1):
                if sampleNode1 is sampleNode2:
                    forward = self.fl / 4.0 > random.random()
                else:
                    forward = self.fl / 2.0 > random.random()
                    if forward == True:
                        c1 = c1.circularize()
                        c2 = c2.circularize()
                gain = not forward
                 
            # case 3) no gain or loss
            else:
                pass

            # do the dcj
            dcjResult = dcj(c1, offset1,
                            c2, offset2,
                            random.randint(0, 1) == 1)

            # do assert checks here
            if gain is True:
                pass
            elif gain is False:
                pass
            else:
                pass
            
            # add the resulting contigs back to the pool
            for res in dcjResult:
                self.pool.insert(res, res.size)

    def __ldEvent(self):
        print "gbg"
        pass

    def __garbageSplitEvent(self):
        pass

    def __garbageSelfEvent(self):
        pass

    def __linearFissionFusionEvent(self):
        pass

    # randomly choose an edge in the garbage (none if we select from
    # the pool instead)
    def __sampleGarbage(self):
        x = random.randint(0, self.pool.size() + self.garbage.size - 1)
        if x < self.garbage.size:
            return x
        return None
            
        
        
    
