#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import copy
from contigSim.src.contig import CircularContig
from contigSim.src.contig import LinearContig

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
        
    def testLinearContigs(self):
        lc1 = LinearContig(10)
        assert lc1.numBases() == 9
        assert lc1.isLinear() == True
        assert lc1.isCircular() == False
        cc1 = lc1.circularize()
        assert cc1.isCircular() == True
        assert cc1.numBases() == 9
        lc1 = cc1.linearize()
        assert lc1.numBases() == 9

        rcut1 = lc1.cutRightOff(7)
        assert rcut1.size == 3
        assert lc1.size == 8

        lcut1 = lc1.cutLeftOff(2)
        assert lc1.size == 6
        assert lcut1.size == 3

        lcut1.joinToLeft(rcut1)
        assert lcut1.size == 5

    def testCircularContigs(self):
        cc1 = CircularContig(5)
        lc1 = cc1.linearize()
        assert cc1.isCircular() == True
        assert lc1.isLinear() == True
        assert lc1.size == 6
        assert cc1.numBases() == 5
        cc2 = cc1.cutOffCircle(0, 3)
        assert cc2.size == 3
        assert cc1.size == 2

        cc1.joinOval(cc2)
        assert cc1.numBases() == 5
        
        
        
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

