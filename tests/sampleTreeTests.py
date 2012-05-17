#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
from contigSim.src.sampleTree import SampleTree

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
        
    def testSampleTreeConstruct(self):
        tree = SampleTree(degree=2)
        assert tree.size() == 0

        tree.insert("cat", 3)
        assert tree.size() == 1
        assert tree.weight() == 3

        tree.insert("bear", 100)
        assert tree.size() == 2
        assert tree.weight() == 103

        tree.insert("rabbit", 100)
        tree.printWeights(tree.root)
        assert tree.size() == 3
        assert tree.weight() == 203

    def testSampleTreeBigger(self):
        tree = SampleTree(degree=5)
        for i in range(0,10000):
            tree.insert(str(i), i)
        assert tree.size() == 10000
        assert tree.weight() == (9999 * 10000) / 2
        

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

