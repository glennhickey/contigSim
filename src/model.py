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
        self.__resetCounts()

    ##################################################################
    # there are five kinds of rates:
    # N: (fixed) number of bases in the model
    # rll: rate for dcj on the bases in the contig pool
    # rld: rate for dcj where one break is in the pool
    #      and the other rate is in the garbage
    # rdd: both in garbage
    # fl: telomere loss modifier
    # fg: telomere gain modifier
    # pgain: dead gain probability
    ##################################################################
    def setParameters(self, N, rll, rld = 0, rdd = 0, fl = 0, fg = 0,
                      pgain = 0):        
        self.eventQueue.reset()
        self.N = N
        self.fl = fl
        self.fg = fg
        self.pgain = pgain

        if rll > 0:
            self.eventQueue.addEventType(N * rll, self.__llEvent)
        if rld > 0:
            self.eventQueue.addEventType(N * rld, self.__ldEvent)
        if rdd > 0:
            self.eventQueue.addEventType(N * rdd, self.__ddEvent)
            
    ##################################################################
    # intitialize the starting state
    # the the contigs will all have the same sizes (modulo rounding)
    # in order to satisfy the input parameters exactly
    ##################################################################
    def setStartingState(self, garbageSize, numLinear, numCircular):
        assert self.N > garbageSize + numLinear + numCircular
        self.pool = SampleTree()

        numGarbage = 0
        if garbageSize > 0:
            garbage = CircularContig(garbageSize)
            garbage.setDead()
            self.pool.insert(garbage, garbage.numBases())
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
                self.pool.insert(contig, contig.numBases())
                added += contig.size
            assert added == linearBases + numLinear
            assert self.pool.size() == numLinear + numGarbage
            assert self.pool.weight() == linearBases + garbageSize

        if numCircular > 0:
            circSize = math.floor(circularBases / numCircular)
            extra = circularBases % numCircular
            added = 0
            for i in range(numCircular):
                size = circSize
                if i < extra:
                    size += 1
                contig = CircularContig(size)
                self.pool.insert(contig, contig.numBases())
                added += contig.size
            assert added == circularBases
            assert self.pool.size() == numLinear + numCircular + numGarbage
            assert self.pool.weight() == circularBases + linearBases + \
            garbageSize

    ##################################################################
    # run the simulation for the specified time
    ##################################################################
    def simulate(self, time):
        self.eventQueue.begin()
        self.__resetCounts()
        while True:
            nextEvent = self.eventQueue.next(time)
            if nextEvent is not None:
                nextEvent()
            else:
                break

    ##################################################################
    # draw (and remove) two random adajcenies and their
    # contigs from the pool (only if they are not dead)
    ##################################################################
    def __drawSamples(self):
        sampleNode1, offset1 = self.pool.uniformSample()
        sampleNode2, offset2 = self.pool.uniformSample()

        # the offset is weighted based on the number of bases
        # we want to translate this into number of edges (splitting)
        # the probability between linear and telomere edges.
        # so for linear contigs with zero offset, we flip a coin to
        # move it to the other side. 
        if sampleNode1.data.isLinear() and offset1 == 0:
            if random.random() < 0.5:
                offset1 = sampleNode1.data.numBases()
        if sampleNode2 is not sampleNode1 and sampleNode2.data.isLinear() and\
           offset2 == 0:
            if random.random() < 0.5:
                offset2 = sampleNode2.data.numBases()

        assert offset1 < sampleNode1.data.size
        assert offset2 < sampleNode2.data.size
        
        return (sampleNode1, offset1, sampleNode2, offset2)

    
    ##################################################################
    #LIVE-LIVE event.  Is normal DCJ operation between two live contigs
    #unless the two breakpoints are identical or on telomeres, in which
    #case fl and fg parameters are used to use fission operations to
    #modifiy the number of telomeres
    ##################################################################
    def __llEvent(self):
        if self.pool.size() == 0 or self.pool.weight() == 1:
            return
        
        # draw (and remove) two random adajcenies and their
        #contigs from the pool (only if they are not dead)
        sampleNode1, offset1, sampleNode2, offset2 = self.__drawSamples()
        c1 = sampleNode1.data
        c2 = sampleNode2.data

        # don't deal with dead contigs in this event
        if c1.isDead() == True or c2.isDead() == True:
            return

        self.pool.remove(sampleNode1)
        if c1 is not c2:
            self.pool.remove(sampleNode2)

        # case 1) gain of telomere
        if sampleNode1 is sampleNode2 and offset1 == offset2:
            return self.__llGain(c1, c2, offset1, offset2)
            
          
        # case 2) loss of telomere
        elif c1.isLinear() and c2.isLinear() and \
                 (offset1 == 0 or offset1 == c1.size - 1) and \
                 (offset2 == 0 or offset2 == c2.size - 1):
            return self.__llLoss(c1, c2, offset1, offset2)

        # case 3) no gain or loss
        self.llCount += 1
        forward = random.randint(0, 1) == 1

        # do the dcj
        dcjResult = dcj(c1, offset1, c2, offset2, forward)
            
        # add the resulting contigs back to the pool
        for res in dcjResult:
            self.pool.insert(res, res.numBases())
            
    ##################################################################
    # Do the fission telomere gain operation (if fg check passes)
    ##################################################################
    def __llGain(self, c1, c2, offset1, offset2):
        # correct "not composite check below"
        if c1.isCircular() or (offset1 != 0 and offset1 != c1.size - 1):
            forward = self.fg > random.random()
            if forward:
                self.fgCount += 1
                dcjResult = dcj(c1, offset1, c2, offset2, forward)
                if c1.isCircular():
                    assert len(dcjResult) == 1 and dcjResult[0].isLinear()
                else:
                    assert len(dcjResult) == 2 and dcjResult[0].isLinear() \
                           and dcjResult[1].isLinear()
                # add the resulting contigs back to the pool
                for res in dcjResult:
                    self.pool.insert(res, res.numBases())
                return

        self.pool.insert(c1, c1.numBases())
        if c2 is not c1:
            self.pool.insert(c2, c2.numBases())
                     
    ##################################################################
    # Do the fission telomer loss operation (if fl check passes)
    ##################################################################
    def __llLoss(self, c1, c2, offset1, offset2):
        if c1 is c2:
            forward = self.fl / 4.0 > random.random()
        else:
            forward = self.fl / 2.0 > random.random()
            if forward == True:
                c1 = c1.circularize()
                c2 = c2.circularize()
        if forward:
            dcjResult = dcj(c1, offset1, c2, offset2, forward)
            self.flCount += 1
            assert c1.isLinear() and c2.isLinear()
            assert len(dcjResult) == 1
            if c1 is not c2:
                assert dcjResult[0].isLinear()
            else:
                assert dcjResult[0].isCircular()
            # add the resulting contigs back to the pool
            for res in dcjResult:
                self.pool.insert(res, res.numBases())
        else:
            self.pool.insert(c1, c1.numBases())
            if c2 is not c1:
                self.pool.insert(c2, c2.numBases())


    ##################################################################
    #LIVE-DEAD (or DEAD-LIVE) event.  One contig is alive and the
    #other is the unique dead contig.  This can result in a loss of
    #live contigs and/or change in number of live bases
    ##################################################################
    def __ldEvent(self):
        if self.pool.size() == 0 or self.pool.weight() == 1:
            return
        
        # draw (and remove) two random adajcenies and their
        #contigs from the pool (only if they are not dead)
        sampleNode1, offset1, sampleNode2, offset2 = self.__drawSamples()
        c1 = sampleNode1.data
        c2 = sampleNode2.data

        # only deal with live / dead contigs in this event
        if (c1.isDead() == c2.isDead()):
            return

        self.pool.remove(sampleNode1)
        if c1 is not c2:
            self.pool.remove(sampleNode2)

        # make sure c1 is alive and c2 is dead
        if c1.isDead():
            c1, c2 = c2, c1
            offset1, offset2 = offset2, offset1

        # do the dcj
        dcjResult = dcj(c1, offset1, c2, offset2, random.randint(0, 1) == 1)

        deadIdx = 0;
        if len(dcjResult) == 2 and \
               random.randint(0, dcjResult[0].size + dcjResult[1].size) >= \
               dcjResult[0].size:
            deadIdx = 1
        dcjResult[deadIdx].setDead(True)

        if len(dcjResult) == 1:
            self.ldLossCount += 1
        else:
            self.ldSwapCount += 1
            
         # add the resulting contigs back to the pool
        for res in dcjResult:
            self.pool.insert(res, res.numBases())
            
    ##################################################################
    #DEAD-DEAD event.  The dead contig rearranges with itself.  pgain
    #is used to decide how oftern this oepration breaks off a new circular
    #live chormosome
    ##################################################################
    def __ddEvent(self):
        if self.pool.size() == 0 or self.pool.weight() == 1:
            return
        
        sampleNode1, offset1, sampleNode2, offset2 = self.__drawSamples()
        c1 = sampleNode1.data
        c2 = sampleNode2.data

        # only deal with dead / dead contigs in this event
        if (c1.isDead() == False or c2.isDead() == False):
            return

        # only support single dead contig
        assert c1 is c2

        # don't know what to do here
        if (offset1 == offset2):
            return        

        self.pool.remove(sampleNode1)
        if c1 is not c2:
            self.pool.remove(sampleNode2)

        #forward means do not cut
        forward = random.random() > self.pgain

        # do the dcj
        dcjResult = dcj(c1, offset1, c2, offset2, forward)

        deadIdx = 0;
        if len(dcjResult) == 2 and \
               random.randint(0, dcjResult[0].size + dcjResult[1].size) \
                            >= dcjResult[0].size:
                    deadIdx = 1
        dcjResult[deadIdx].setDead(True)

        if forward:
            self.ddSwapCount += 1
            assert len(dcjResult) == 1
        else:
            self.ddGainCount += 1
            assert len(dcjResult) == 2
        
         # add the resulting contigs back to the pool
        for res in dcjResult:
            self.pool.insert(res, res.numBases())

          
    ##################################################################
    # all counters set to zero.  
    ##################################################################
    def __resetCounts(self):
        self.llCount = 0
        self.fgCount = 0
        self.flCount = 0
        self.ldLossCount = 0
        self.ldSwapCount = 0
        self.ddGainCount = 0
        self.ddSwapCount = 0

        
        
    
