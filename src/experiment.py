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

def main(argv=None):
    if argv is None:
        argv = sys.argv

    exp = Experiment()
    t = 1000000
    N = 10000
    exp.addParameterSet(t, N, 0.1 / N)
    exp.addStartingState(10, 10, 10)
    exp.run(1, 1)
    print exp.results
    

    
    
if __name__ == "__main__":
    sys.exit(main())
