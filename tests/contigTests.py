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
        assert cc1.numBases() == 11
        lc1 = cc1.linearize()
        assert lc1.numBases() == 9

        l,r = lc1.cut(7)
        assert l.size == 7
        assert r.size == 2
        lc1 = l.joinToLeft(r)
        assert lc1.size == 10

    def testCircularContigs(self):
        cc1 = CircularContig(5)
        lc1 = cc1.linearize()
        assert cc1.isCircular() == True
        assert lc1.isLinear() == True
        assert lc1.size == 4
        assert cc1.numBases() == 5

        l,r  = cc1.cut(0, 3)
        assert l.size == 2
        assert r.size == 3

        cc1 = l.join(r)
        assert cc1.numBases() == 5
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

