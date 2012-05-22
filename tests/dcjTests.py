#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import copy

from contigSim.src.eventQueue import EventQueue
from contigSim.src.contig import CircularContig
from contigSim.src.contig import LinearContig
from contigSim.src.dcj import dcj

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
        
    def testDcjLinear(self):
        cont = LinearContig(10)
        res = dcj(cont, 2, cont, 6, True)
        assert len(res) == 1
        assert res[0].size == 10
        assert res[0].isLinear() == True

        res = dcj(cont, 2, cont, 6, False)
        assert len(res) == 2
        assert res[0].isLinear() == True
        assert res[1].isCircular() == True
        assert res[0].size == 6
        assert res[1].size == 4

        res = dcj(cont, 1, cont, 0, False)
        assert len(res) == 2
        assert res[0].isLinear() == True
        assert res[1].isCircular() == True
        assert res[0].size == 9
        assert res[1].size == 1

        res = dcj(cont, 1, cont, 1, False)
        assert len(res) == 1
        assert res[0].isLinear()
        assert res[0].size == cont.size

        res = dcj(cont, 1, cont, 1, True)
        assert len(res) == 2
        assert res[0].isLinear() and res[1].isLinear()
        assert res[0].size + res[1].size == cont.size + 1

    def testDcjLinearLinear(self):
        c1 = LinearContig(100)
        c2 = LinearContig(50)
        res = dcj(c1, 30, c2, 20, True)
        assert len(res) == 2
        assert res[0].size == 51
        assert res[1].size == 99

        res = dcj(c1, 30, c2, 20, False)
        assert len(res) == 2
        assert res[0].size + res[1].size == 150
        
    def testDcjLinearCircular(self):
        c1 = LinearContig(10)
        c2 = CircularContig(5)
        res = dcj(c1, 3, c2, 0, True)
        assert len(res) == 1
        assert res[0].size == 15
        assert res[0].isLinear() == True

        res = dcj(c1, 3, c2, 3, False)
        assert len(res) == 1
        assert res[0].size == 15
        assert res[0].isLinear() == True

    def testDcjCircular(self):
        cont = CircularContig(10)
        res = dcj(cont, 3, cont, 4, True)
        assert len(res) == 1
        assert res[0].isCircular() == True
        assert res[0].size == 10

        res = dcj(cont, 2, cont, 9, False)
        assert len(res) == 2
        assert res[0].isCircular() == True
        assert res[0].size == 7
        assert res[1].isCircular() == True
        assert res[1].size == 3

        res = dcj(cont, 2, cont, 2, False)
        assert len(res) == 1
        assert res[0].size == cont.size
        assert type(res[0]) == type(cont)

        res = dcj(cont, 2, cont, 2, True)
        assert len(res) == 1
        assert res[0].size == cont.size
        assert res[0].isCircular() == False


    def testDcjCircularCircular(self):
        c1 = CircularContig(33)
        c2 = CircularContig(15)
        res = dcj(c1, 10, c2, 11, True)
        assert len(res) == 1
        assert res[0].isCircular() == True
        assert res[0].size == 48

        res = dcj(c1, 1, c2, 11, False)
        assert len(res) == 1
        assert res[0].isCircular() == True
        assert res[0].size == 48

    def testDcjCircularLinear(self):
        c2 = LinearContig(10)
        c1 = CircularContig(5)
        res = dcj(c1, 3, c2, 0, True)
        assert len(res) == 1
        assert res[0].size == 15
        assert res[0].isLinear() == True

        res = dcj(c1, 3, c2, 3, False)
        assert len(res) == 1
        assert res[0].size == 15
        assert res[0].isLinear() == True
    
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

