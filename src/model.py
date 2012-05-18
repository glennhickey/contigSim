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
        self.garbage = None

    # there are five kinds of rates:
    # poolRate: rate for dcj on the bases in the contig pool
    # poolGarbageRate: rate for dcj where one break is in the pool
    #      and the other rate is in the garbage
    # garbageSplitRate: rate at which the garbage cuts off a piece of itself
    #      to form a new contig (which gets added to pool)
    # garbageSelfRate: rearrangements which just shuffle the garbage
    # linearFissionFusionRate: split or merge linear contigs
    def setRates(self, poolRate, poolGarbageRate = 0, garbageSplitRate = 0,
                 garbageSelfRate = 0, linearFissionFusionRate = 0):
        self.eventQueue.reset()
        if poolRate > 0:
            self.eventQueue.addEventType(poolRate,
                                         self.__poolEvent)
        if poolGarbageRate > 0:
            self.eventQueue.addEventType(poolGarabgeRate,
                                         self.__poolGarbageEvent)
        if garbageSplitRate > 0:
            self.eventQueue.addEventType(poolGarabgeRate,
                                         self.__garbageSplitEvent)
        if garbageSelfRate > 0:
            self.eventQueue.addEventType(poolGarabgeRate,
                                         self.__garbageSelfEvent)
        if linearFissionFusionRate > 0:
            self.eventQueue.addEventType(poolGarabgeRate,
                                         self.__linearFissionFusionEvent)

    # intitialize the starting state
    def setStartingState(self, nBases, garbageSize, numLinear, numCircular):
        assert nBases > garbageSize + numLinear + numCircular
        if garbageSize > 0:
            self.garbage = CircularContig(garbageSize)
        else:
            self.garbage = None
        linearBases = math.floor((nBases - garbageSize) / 2.0)
        circularBases = math.ceil((nBases - garbageSize) / 2.0)

        linSize = math.floor(linearBases / numLinear)
        extra = linearBases % numLinear
        added = 0
        for i in range(numLinear):
            size = linSize
            if i < extra:
                size += 1
            # plus 1 because the number of adjacencies is 1 + number of bases
            contig = LinearContig(size + 1)
            self.pool.insert(contig, contig.size)
            added += contig.size
        assert added == linearBases + numLinear
        assert self.pool.size() == numLinear
        assert self.pool.weight() == linearBases + numLinear

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
        assert self.pool.size() == numLinear + numCircular
        assert self.pool.weight() == circularBases + linearBases + numLinear

    # run the simulation for the specified time
    def simulate(self, time):
        self.eventQueue.begin()
        while True:
            nextEvent = self.eventQueue.next(time)
            if nextEvent is not None:
                nextEvent()
            else:
                break
    
    def __poolEvent(self):
        print "pool"
        pass

    def __poolGarbageEvent(self):
        print "gbg"
        pass

    def __garbageSplitEvent(self):
        pass

    def __garbageSelfEvent(self):
        pass

    def __linearFissionFusionEvent(self):
        pass

            
        
        
    
