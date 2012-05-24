#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import copy
import random
from collections import defaultdict

""" In order to quickly sample contigs (uniformly based on their weights)
we keep them in a b-tree.  This way sampling can be done in logN, as can
updates due to rearrangements.  

"""

class SampleTreeNode(object):
    def __init__(self, parent):
        self.weight = 0
        self.count = 0
        self.parent = parent
        self.children = []
        self.data = None

class SampleTree(object):
    def __init__(self, degree=4):
        self.degree = degree
        self.root = SampleTreeNode(None)
        assert self.degree > 1

    # find a free slot in the tree
    def __findSlot(self, node):
        # case 1: internal node has room for another child
        if node.data is None and len(node.children) < self.degree:
            return node
        # case 2: on leaf: turn into an internal node and push leaf down
        elif node.data is not None:
            assert len(node.children) == 0
            child = SampleTreeNode(node)
            child.data = node.data
            child.weight = node.weight
            child.count = node.count
            node.data = None
            node.children.append(child)
            return node
        # case 3: recurse down full internal node
        else:
            child = min(node.children, key=lambda x: x.count)
            return self.__findSlot(child)

    # update weight and count values on parents of node
    def __updateUpwards(self, node):
        if node is None:
            return
        assert node.data is None
        node.weight = reduce(lambda x,y: x + y.weight, node.children, 0)
        node.count = reduce(lambda x,y: x + y.count, node.children, 0)
        self.__updateUpwards(node.parent)

    # insert a new leaf node with given data and weight
    def insert(self, data, weight):
        parent = self.__findSlot(self.root)
        assert len(parent.children) < self.degree and parent.data is None
        newNode = SampleTreeNode(parent)
        newNode.data = data
        newNode.weight = weight
        newNode.count = 1
        parent.children.append(newNode)
        self.__updateUpwards(newNode.parent)

    # remove a given leaf node
    def remove(self, node):
        assert len(node.children) == 0 and node.data is not None
        assert node.parent is not None
        node.parent.children = [x for x in node.parent.children if x != node]
        self.__updateUpwards(node.parent)    

    # how many data elememnts are in the tree
    def size(self):
        return self.root.count

    # what is the total weight of the data elements in the tree
    # the probability of selecting an element is its weight over the total
    def weight(self):
        return self.root.weight

    # uniformly sample a data element based on its weight
    def uniformSample(self, node=None):
        if node is None:
            node = self.root
        if node.weight == 0:
            return None
        x = random.randint(0, node.weight - 1)
        tally = int(0)
        for child in node.children:
            if x < tally + child.weight:
                if len(child.children) == 0:
                    assert x - tally < child.weight
                    return (child, x - tally)
                else:
                    return self.uniformSample(child)
            tally += child.weight
        assert False

    # iterate through the nodes containing data elements
    #(stored in leaves) in the tree
    def nodes(self, node=None):
        if node is None:
            node = self.root
        yield node
        for child in node.children:
            for recChild in self.nodes(child):
                yield recChild
                
    # iterate through the data elements (stored in leaves) in the tree
    def dataElements(self, node=None):
        if node is None:
            node = self.root
        for i in self.nodes(node):
            if i.data is not None:
                yield i.data
        
    # git a histogram of the node weights of data elements with whose
    # types are instances of the given dataType
    def histogram(self, binSize = 1, dataType=None, checkFn = None):
        hist = defaultdict(int)
        for node in self.nodes():
            if node.data is not None and\
               (dataType is None or issubclass(type(node.data), dataType)) and\
               (checkFn is None or checkFn(node.data) == True):
                bin = int(node.weight) / int(binSize)
                hist[bin] += 1
        return hist

    def printWeights(self, node):
        print "%d " % node.weight
        for child in node.children:
            self.printWeights(child)
        
    
    
    
