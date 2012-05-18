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
        res = dcj(cont, 2, 6, True)
        assert len(res) == 1
        assert res[0].size == 10
        assert res[0].isLinear() == True

        res = dcj(cont, 2, 6, False)
        assert len(res) == 2
        assert res[0].isLinear() == True
        assert res[1].isCircular() == True
        assert res[0].size == 6
        assert res[1].size == 4

        res = dcj(cont, 1, 0, False)
        assert len(res) == 2
        assert res[0].isLinear() == True
        assert res[1].isCircular() == True
        assert res[0].size == 9
        assert res[1].size == 1

    def testDcjLinearLinear(self):
        c1 = LinearContig(100)
        c2 = LinearContig(50)
        res = dcj(c1, 30, 20, True, c2)
        assert len(res) == 2
        assert res[0].size == 51
        assert res[1].size == 99
        
    def testDcjLinearCircular(self):
        
        pass
        
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

