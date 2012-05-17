#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txtimport unittest

import unittest
import sys
import os
from contigSim.tests.sampleTreeTests import TestCase as sampleTreeTest
from contigSim.tests.eventQueueTests import TestCase as eventQueueTest

def allSuites(): 
    allTests =unittest.TestSuite(
        (unittest.makeSuite(sampleTreeTest, 'test'),
         unittest.makeSuite(eventQueueTest, 'test')))
    return allTests
        
def main():    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
                
