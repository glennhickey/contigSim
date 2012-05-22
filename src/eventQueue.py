#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import copy
import random
import math
from heapq import heappush, heappop

""" Generate events over time.  There can be different types of events
with different (exponential) rates.  

"""
class EventQueue(object):
    def __init__(self):
        self.time = 0
        self.rates = dict()
        self.heap = []

    # zap everything
    def reset(self):
        self.time = 0
        self.rates = dict()
        self.heap  = []

    # add a new event (with an exponential rate and a unique name)
    def addEventType(self, rate, name):
        assert name not in self.rates
        self.rates[name] = rate

    # begin the simulation at time=0
    def begin(self):
        self.time = 0
        self.heap = []
        for item in self.rates.items():
            delta = random.expovariate(item[1])
            heappush(self.heap, (delta, item[0]))

    # move clock forward to next event and return its name
    def next(self, maxTime=sys.maxint):
        if len(self.heap) == 0:
            return None
        item = heappop(self.heap)
        name = item[1]
        assert self.time <= item[0]
        self.time = item[0]
        if self.time > maxTime:
            self.time = maxTime
            return None
        else:
            nextTime = random.expovariate(self.rates[name]) + self.time
            heappush(self.heap, (nextTime, name))
            return name
