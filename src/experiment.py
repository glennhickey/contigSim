#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from model import Model
from sampleTree import SampleTree


# framework for generating experimental results from the simulation,
# over differen sets of parameters
class Experiment(object):
    def __init__(self):
        self.replicates = 1
        self.parameterSpace = []
        self.startingStateSpace = []
        self.results = dict()
        self.binSize = 1

    def addParameterSet(self, t, N, rll, rld = 0, rdd = 0, fl = 0, fg = 0,
                        pgain = 0):
        self.parameterSpace.append((t, N, rll, rld, rdd, fl, fg, pgain))

    def addStartingState(self, garbageSize, numLinear, numCircular):
        self.startingStateSpace.append((garbageSize, numLinear, numCircular))

    def setNumReplicates(self, n):
        self.replicates = n

    def reset():
        self.parameterSpace = []
        self.startingStateSpace = []
        self.results = dcit
        self.replicates = 1

    def run(self, replicates = 1, binSize = 1):
        self.replicates = replicates
        self.binSize = binSize
        for params in self.parameterSpace:
            for state in self.startingStateSpace:
                self.__runInstance(params, state)

    def __extractResultsFromModel(self, model):
        res = dict()
        res["overall"] = model.pool.histogram(binSize=self.binSize)
        res["dead"] = model.pool.histogram(binSize=self.binSize,
            checkFn = lambda x : x.isDead() == True)
        res["alive"] =  model.pool.histogram(binSize=self.binSize,
            checkFn = lambda x : x.isDead() == False)
        res["aliveLinear"] = model.pool.histogram(binSize=self.binSize,
            checkFn = lambda x : not x.isDead() and x.isLinear())
        res["aliveCircular"] = model.pool.histogram(binSize=self.binSize,
            checkFn = lambda x : not x.isDead() and x.isCircular())
        return res

    def __runInstance(self, parameters, startState):
        for rep in xrange(0, self.replicates):
            model = Model()
            model.setParameters(parameters[1], parameters[2], parameters[3],
                                parameters[4], parameters[5], parameters[6],
                                parameters[7])
            model.setStartingState(startState[0], startState[1], startState[2])
            model.simulate(parameters[0])
            results = self.__extractResultsFromModel(model)
            key = parameters + startState
            if rep is 0:
                self.results[key] = []
            self.results[key].append(results)
            print (model.llCount,
                   model.fgCount,
                   model.flCount,
                   model.ldLossCount,
                   model.ldSwapCount,
                   model.ddGainCount,
                   model.ddSwapCount)

def avgHistogram(results, cat, N):
    table = defaultdict(int)
    for rep in results:
        res = rep[cat]
        for key,value in res.items():
            table[key] += value
    for key,value in table.items():
        table[key] = float(table[key]) / float(len(results))
    return table

def cumulative(table):
    xAxis = []
    yAxis = []    
    for key,value in table.items():
            xAxis.append(key)
            yAxis.append(value * 1)
    return (xAxis, yAxis)
    
def doPlot(ctable, ltable, title):
    #n, bins, patches = plt.hist(table, 50, normed=1, facecolor='green',
    #                            alpha=0.75)

    # add a 'best fit' line
    #y = mlab.normpdf( bins, mu, sigma)
    #l = plt.plot(bins, y, 'r--', linewidth=1)
    xAxis, yAxis = cumulative(ctable)
    xAxis2, yAxis2 = cumulative(ltable)
    plt.plot(xAxis, yAxis, 'r+', linewidth=1,markersize=25)
    plt.plot(xAxis2, yAxis2, 'g+', linewidth=1,markersize=10)
    #plt.hist(xrange(len(table)), table)
    plt.xlabel('i')
    plt.ylabel('iEi')
    plt.title(title)
#    plt.axis([0, 100000, 0, 100000])
    plt.grid(True)
    print xAxis, yAxis
    
    plt.show()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    exp = Experiment()
    t = 10000
    N = 3000000000
    binSize = 1000000
    replicates = 50
    exp.addParameterSet(t, N, 1.0 / N, 0, fl = 0.00, fg = 0.00, pgain = 0.00)
    exp.addStartingState(0, 25, 0)
    exp.run(replicates, binSize)

    # there iterate for each combination of starting state and parameters
    # note there is only one for now.
    # each result is a name - table pair
    for result in exp.results.items():

        # make unique filename as function of parameters
        fname = str(result[0])
        fname = fname.replace(" ", "").replace("(", "").replace(")", "")
        fname = fname.replace(",", "_")
        # add extensions (pdf?)
        fname += ".pdf"

        # the results tables are lists of replicates
        # use the crappy avgHistogram method to compute
        # the avearge (mean) of each table
        ctable = avgHistogram(result[1], "aliveCircular", N)
        ltable = avgHistogram(result[1], "aliveLinear", N)
        dtable = avgHistogram(result[1], "dead", N)

        # run through some basic counts for debugging purposes
        numLinearBases = 0
        numLinearContigs = 0
        for key,value in ltable.items():
            numLinearBases += key * value
            numLinearContigs += value
        numCircularBases = 0
        numCircularContigs = 0
        for key,value in ctable.items():
            numCircularBases += key * value
            numCircularContigs += value
        numDeadBases = 0
        numDeadContigs = 0
        for key,value in dtable.items():
            numDeadBases += key * value
            numDeadContigs += value

        print "linear: contigs=%d bases=%d" % (numLinearContigs,
                                               numLinearBases)
        print "circular: contigs=%d bases=%d" % (numCircularContigs,
                                                 numCircularBases)

        print "dead: contigs=%d bases=%d" % (numDeadContigs,
                                             numDeadBases)

        # use current travesty to plot just the circlular and linear
        doPlot(ctable, ltable, fname)

if __name__ == "__main__":
    sys.exit(main())
