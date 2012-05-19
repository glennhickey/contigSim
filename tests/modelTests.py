#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import copy

from contigSim.src.model import Model
from contigSim.src.contig import LinearContig
from contigSim.src.contig import CircularContig

from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        
    def testModelInit(self):
        model = Model()
        assert model.pool.size() == 0

        model.setRates(0.1)
        model.setStartingState(100, 0, 21, 3)
        assert model.pool.size() == 24
        assert model.pool.weight() == 121

        linCount = 0
        cirCount = 0
        for contig in model.pool.dataElements():
            if type(contig) == LinearContig:
                linCount += 1
            elif type(contig) == CircularContig:
                cirCount += 1
            else:
                assert False

        assert linCount == 21
        assert cirCount == 3

        model.setStartingState(100, 1, 3, 0)
        assert model.pool.size() == 3
        assert model.pool.weight() == 102

        model.setStartingState(55, 11, 0, 4)
        assert model.pool.size() == 4
        assert model.pool.weight() == 44

    def testSimulateCircularPool(self):
        model = Model()
        model.setRates(0.1)
        model.setStartingState(100, 0, 0, 2)
        model.simulate(100)
        assert model.eventQueue.time == 100
        
    
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

