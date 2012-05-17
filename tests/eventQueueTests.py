#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
from contigSim.src.eventQueue import EventQueue

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
        
    def testEventQueueBuild(self):
        eq = EventQueue()
        eq.addEventType(0.01, "substitution")
        eq.addEventType(0.001, "duplication")
        eq.addEventType(0.0005, "fusion")

        eq.begin()
        events = []
        counts = dict()
        counts["substitution"] = 0
        counts["duplication"] = 0
        counts["fusion"] = 0
        
        for i in range(0, 10000):
            eventName = eq.next()
            counts[eventName] += 1

        assert counts["substitution"] >= counts["duplication"]
        assert  counts["duplication"] >= counts["fusion"]

        print "time %f" % eq.time 
        print "%d subs;  %d dups;  %d fus" % (counts["substitution"],
                                          counts["duplication"],
                                           counts["fusion"])
        
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

