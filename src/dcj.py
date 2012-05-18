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

from contig import CircularContig
from contig import LinearContig
from contig import Contig

# general interface to a djc operation (of which there are six possible kinds)
# the positions are the target edges to be cut, and forward=True specifes
# they are reattached left-to-left (False would be left-to-right)
def dcj(cont1, pos1, pos2, forward=True, cont2=None):
    t1 = type(cont1)
    t2 = type(cont2)
    assert issubclass(t1, Contig)
    assert cont2 is None or issubclass(t2, Contig)

    if t1 == LinearContig:
        if cont2 is None:
            return __dcj_linear(cont1, pos1, pos2, forward)
        elif t2 == LinearContig:
            return __dcj_linear_linear(cont1, pos1, pos2, forward, cont2)
        else:
            return __dcj_linear_circular(cont1, pos1, pos2, forward, cont2)

    if t1 == CircularContig:
        if cont2 is None:
            return __dcj_circular(cont1, pos1, pos2, forward)
        elif t2 == LinearContg:
            return __dcj_circular_circular(con1, pos1, pos2, forward, cont2)
        else:
            return __dcj_circular_linear(con1, pos1, pos2, forward, cont2)

    assert False

# dcj on a single linear contig
# case 1) both breaks are on same edge: (ERROR)
# case 2) join in "forward sense" returns A-BC
# case 3) join in "reverse sense" returns AC, circle(B)"
def __dcj_linear(cont, pos1, pos2, forward):
    assert pos1 != pos2
    p1 = min(pos1, pos2)
    p2 = max(pos1, pos2)
    left, temp = cont.cut(p1)
    middle, right = temp.cut(p2 - left.size - 1)
    if forward:
        temp = left.joinToRight(middle, forward=False)
        return (temp.joinToRight(right, forward=True),)
    else:
        linCont = left.joinToRight(right, forward=True)
        cirCont = middle.circularize()
        return (linCont, cirCont)

def __dcj_linear_linear(cont1, pos1, pos2, forward, cont2):
    pass

def __dcj_linear_circular(cont1, pos1, pos2, forward, cont2):
    pass

def __dcj_circular(cont, pos1, pos2, forward):
    pass

def __dcj_circular_circular(cont1, pos1, pos2, forward, cont2):
    pass

def __dcj_circular_linear(cont1, pos1, pos2, forward, cont2):
    pass
