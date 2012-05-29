#!/usr/bin/env python
"""
runSim.py
28 May 2012
Glenn Hickey, hickey (a) soe ucsc edu
Dent Earl, dearl (a) soe ucsc edu


"""
#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
import argparse
import cPickle
import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import matplotlib.mlab as mlab

from contigSim.src.model import Model
from contigSim.src.sampleTree import SampleTree
from contigSim.src.experiment import Experiment

def initOptions():
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('--savePlot', type=str, help='Location to save pdf plot.')
    parser.add_argument('--replicates', type=int, default=50, 
                        help='Number of replicates to run. default=%(default)s')
    parser.add_argument('--binSize', type=int, default=1000000,
                        help='Size of bins. default=%(default)s')
    parser.add_argument('--N', type=int, default=3000000000,
                        help='Number of bases. default=%(default)s')
    parser.add_argument('--t', type=float, default=10000,
                        help='Time. default=%(default)s')
    parser.add_argument('--linearY', default=False, action='store_true',
                        help='Plot y-axis in linear scale. default=%(default)s')
    parser.add_argument('--countY', default=False, action='store_true',
                        help='Y-axis shows counts per bin instead of frequency per bin. default=%(default)s')
    parser.add_argument('--saveSim', type=str, help='Location to save pickle.')
    parser.add_argument('--loadSim', type=str, help='Location to load pickle.')
    parser.add_argument('--alpha', type=float, default=0.3, 
                        help='Alpha transparency for plot, [0, 1]. default=%(default)s')
    return parser
def checkOptions(args, parser):
    if args.saveSim is not None and args.loadSim is not None:
        parser.error('Error, you cannot invoke both --loadSim and --saveSim. Pick one.')
    if args.loadSim is not None:
        if not os.path.exists(args.loadSim):
            parser.error('Error, file --loadSim %s does not exist!' % args.loadSim)
def packData(obj, filename):
    """packData() stores an object in filename
    """
    f = open(filename, 'wb')
    cPickle.dump(obj, f, 2) # 2 is the format protocol, 2 = binary     
    f.close()
def unpackData(file):
   """unpackData() opens up the pickle of the last run and pulls out
   all the relevant data.                                           
   """
   f = open(file, 'r')
   obj = cPickle.load(f)
   f.close()
   return obj
def avgHistogram(results, cat, args):
    table = defaultdict(int)
    for rep in results:
        res = rep[cat]
        for key,value in res.items():
            table[key] += value
    for key,value in table.items():
        if args.countY:
            table[key] = float(table[key])
        else:
            table[key] = float(table[key]) / float(len(results))
    return table

def cumulative(table):
    xAxis = []
    yAxis = []    
    for key,value in table.items():
            xAxis.append(key)
            yAxis.append(value * 1)
    return (xAxis, yAxis)
def initImage(width, height, filename, args):
    """
    initImage takes a width and height and returns both a fig and pdf object.
    """
    pdf = pltBack.PdfPages(filename)
    fig = plt.figure(figsize=(width, height), dpi=300, facecolor='w')
    return (fig, pdf)
def initAxis(fig, args):
    args.axLeft = 0.1
    args.axWidth = 0.85
    args.axBottom = 0.2
    args.axHeight = 0.7
    ax = fig.add_axes([args.axLeft, args.axBottom,
                       args.axWidth, args.axHeight])
    return ax
def writeImage(fig, pdf, args):
    fig.savefig(pdf, format = 'pdf')
    pdf.close()
def doPlot(ctable, ltable, dctable, dltable, title, args):
    fig, pdf = initImage(9.0, 4.0, title, args)
    ax = initAxis(fig, args)
    drawData(ax, ctable, ltable, dctable, dltable, title, args)
    if args.savePlot is None:
        plt.show()
    else:
        writeImage(fig, pdf, args)
def extractPlottables(table, binSize):
    keylist = table.keys()
    keylist.sort()
    x = []
    y = []
    for key in keylist:
        x.append(key)
        y.append(table[key])
    # center the x values in the middle of their bin
    x = (np.array(x) * binSize) - binSize / 2.0
    y = np.array(y)
    return x, y 
def drawData(ax, ctable, ltable, dctable, dltable, title, args):
    cx, cy = extractPlottables(ctable, args.binSize)
    lx, ly = extractPlottables(ltable, args.binSize)
    dcx, dcy = extractPlottables(dctable, args.binSize)
    dlx, dly = extractPlottables(dltable, args.binSize)
    colorList = ['#1f77b4', # dark blue
                 '#aec7e8', # light blue
                 '#ff7f0e', # bright orange
                 '#ffbb78', # light orange
                 '#4B4C5E', # dark slate gray
                 '#9edae5', # light blue 
                 '#7F80AB', # purple-ish slate blue
                 '#c7c7c7', # light gray
                 '#9467bd', # dark purple
                 '#c5b0d5', # light purple
                 '#d62728', # dark red
                 '#ff9896', # light red
                 ]
    plotlist = []
    # there is a bug in our version fo matplotlib that wont allow us to set markeredgecolor='none'
    plotlist.append(plt.plot(cx, cy, color=colorList[0], linestyle='none', marker='.', 
                             markeredgecolor=colorList[0], markeredgewidth=0, linewidth=0.5,
                             markersize=10.0, alpha=args.alpha)[0])
    plotlist.append(plt.plot(lx, ly, color=colorList[2], linestyle='none', marker='.', 
                             markeredgecolor=colorList[2], markeredgewidth=0, linewidth=0.5,
                             markersize=10.0, alpha=args.alpha)[0])
    plotlist.append(plt.plot(dcx, dcy, color=colorList[10], linestyle='none', marker='.', 
                             markeredgecolor=colorList[10], markeredgewidth=0, linewidth=0.5,
                             markersize=10.0, alpha=args.alpha)[0])
    plotlist.append(plt.plot(dlx, dly, color=colorList[9], linestyle='none', marker='.', 
                             markeredgecolor=colorList[9], markeredgewidth=0, linewidth=0.5,
                             markersize=10.0, alpha=args.alpha)[0])

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    if xmin < -1.0:
        plt.xlim(0 - .02 * xmax, xmax * 1.02)
    else:
        plt.xlim(xmax - ((xmax - xmin) * 1.02), xmax * 1.02)
    plt.title(title)
    for loc, spine in ax.spines.iteritems():
        if loc in ['left','bottom']:
            # outward by 10 points
            spine.set_position(('outward', 10)) 
        elif loc in ['right','top']:
            # don't draw spine               
            spine.set_color('none') 
        else:
            raise ValueError('unknown spine location: %s' % loc)
    # turn off ticks where there is no spine
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if not args.linearY:
        ax.set_yscale('log')
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(1, 10)))
    plt.xlabel('Binned Contig Size')
    if args.countY:
        plt.ylabel('Count')
    else:
        plt.ylabel('Frequency per Bin')
    leg = plt.legend(plotlist, ['Circular Contigs', 'Linear Contigs', 'Dead Circular Contigs', 'Dead Linear Contigs'],
                     loc=1, numpoints=1)
    plt.setp(leg.get_texts(), fontsize='x-small') # legend fontsize
    leg._drawFrame = False

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = initOptions()
    args = parser.parse_args()
    checkOptions(args, parser)

    if args.loadSim is None:
        exp = Experiment()
        exp.addParameterSet(args.t, args.N, rll=1.0 / args.N,
                            rld=0, rdd=0,
                            fl = 0, fg = 0, pgain = 0.00)        
        exp.addParameterSet(args.t, args.N, rll=1.0 / args.N,
                            rld= 0.1/ args.N, rdd= 0.1 / args.N,
                            fl = 0.5, fg = 0.5, pgain = 0.00)
        exp.addParameterSet(args.t, args.N, rll=1.0 / args.N,
                            rld= 0.1/ args.N, rdd= 0.1 / args.N,
                            fl = 0.5, fg = 0.5, pgain = 0.5)
        exp.addStartingState(0, 25, 0)
        exp.addStartingState(0, 0, 25)
        if args.N > 3000025:
            exp.addStartingState(3000000, 25, 0)
            exp.addStartingState(3000000, 0, 25)
            exp.addStartingState(3000000, 10, 10)
        
        exp.run(args.replicates, 1)
    else:
        exp = unpackData(args.loadSim)
    if args.saveSim is not None:
        packData(exp, args.saveSim)

    # there iterate for each combination of starting state and parameters
    # note there is only one for now.
    # each result is a name - table pair
    for result in exp.results.items():
        # make unique filename as function of parameters
        fname = "t%d_N%d_rll%.2f_rld%.2f_rdd_%.2f_fl%.2f_fg%.2f_pgain_%.2f__gbg%d_nl%d_nc%d" % result[0]
        # add extensions (pdf?)
        txtname = fname + ".txt"
        fname += ".pdf"

        # the results tables are lists of replicates
        # use the crappy avgHistogram method to compute
        # the avearge (mean) of each table
        ctable = avgHistogram(result[1], "aliveCircular", args)
        ltable = avgHistogram(result[1], "aliveLinear", args)
        dctable = avgHistogram(result[1], "deadCircular", args)
        dltable = avgHistogram(result[1], "deadLinear", args)

        # run through some basic counts for debugging purposes
        numLinearBases = 0
        numLinearContigs = 0
        for key,value in ltable.items():
            numLinearBases += key * value
            numLinearContigs += value
        numCircularBases = 0
        numCircularContigs = 0
        for key,value in ctable.items():
            numCircularBases += key * value
            numCircularContigs += value
        numDeadLinearBases = 0
        numDeadLinearContigs = 0
        for key,value in dltable.items():
            numDeadLinearBases += key * value
            numDeadLinearContigs += value
        numDeadCircularBases = 0
        numDeadCircularContigs = 0
        for key,value in dctable.items():
            numDeadCircularBases += key * value
            numDeadCircularContigs += value

        print "linear: contigs=%d bases=%d" % (numLinearContigs,
                                               numLinearBases)
        print "circular: contigs=%d bases=%d" % (numCircularContigs,
                                                 numCircularBases)
        print "deadLinear: contigs=%d bases=%d" % (numDeadLinearContigs,
                                                   numDeadLinearBases)
        print "deadCircular: contigs=%d bases=%d" % (numDeadCircularContigs,
                                                     numDeadCircularBases)

        log = open(txtname, 'w')
        log.write("linear: contigs=%d bases=%d\n" % (numLinearContigs,
                                                     numLinearBases))
        log.write("circular: contigs=%d bases=%d\n" % (numCircularContigs,
                                                       numCircularBases))
        log.write("deadLinear: contigs=%d bases=%d" % (numDeadLinearContigs,
                                                       numDeadLinearBases))
        log.write("deadCircular: contigs=%d bases=%d" % (numDeadCircularContigs,
                                                         numDeadCircularBases))
        log.close()

        # use dent's awesome functinos to plot the results
        doPlot(ctable, ltable, dctable, dltable, fname, args)

if __name__ == "__main__":
    sys.exit(main())
