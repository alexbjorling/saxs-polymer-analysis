#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import Chains
import time
import sys

# default parameters
number = 1
length = 50
box = 10
maxAngle = 90.0
beta = 0.0
outputFile = 'out'
outputFreq = 10
nSteps = 1000

# usage
def usage():
    b = '\n                        '
    print "\nUsage:\n"
    print "chainSimulation.py"
    print "   -surface:         if present, simulates the chains grafted on a hard but" + b + "inert surface (default no)"
    print "   -number <n>:      simulates a collection of <n> chains (default 1)"
    print "   -length <n>:      makes coils of <n> beads, with contour lengths <n>-1" + b + "(default 50)"
    print "   -box <d>:         grafts the chains on a square of side length <d>" + b + "(default 10)"
    print "   -maxAngle <a>:    the maximum bend angle at each bead in degrees where" + b + "<a>=0.0 is a straight rod (default 90.0)"
    print "   -beta <b>:        the effective bond strength <b> = energy / (kT), where" + b + "the well width is 0.2 (default 0.0)"
    print "   -outputFile <s>:  the base name of the outputfile (default: 'out')"
    print "   -outputFreq <n>:  produces output every <n> simulation steps (default: 10)"
    print "   -steps <n>:       runs a simulation of <n> steps, including steps where the" + b + "change is rejected (default 1000)"
    print ""
    exit()

# parse arguments
if '-help' in sys.argv: usage()
surface = '-surface' in sys.argv
def parse(argv, key, default):
    if key in argv:
        return argv[argv.index(key) + 1]
    else:
        return default
number     = int(parse(sys.argv, '-number', number))
length     = int(parse(sys.argv, '-length', length))
box        = float(parse(sys.argv, '-box', box))
maxAngle   = float(parse(sys.argv, '-maxAngle', maxAngle)) * np.pi / 180.0
beta       = float(parse(sys.argv, '-beta', beta))
outputFile = parse(sys.argv, '-outputFile', outputFile)
outputFreq = int(parse(sys.argv, '-outputFreq', outputFreq))
nSteps     = int(parse(sys.argv, '-steps', nSteps))
if (number > 1) and not surface:
    print "Simulating more than one chain without a surface makes no sense!\n"
    exit()

# start the timer and initialize the Chains instance
t0 = time.time()
chains = Chains.Chains(number=number, length=length, box=box, maxAngle=maxAngle, beta=beta, surface=surface, outFile=outputFile+'.pdb')

# set up the output file
try: os.remove(outputFile + '.pdb')
except: pass
chains.dump(append=True)

# do the simulation
goodSteps, badSteps, oldBonds = 0, 0, 0
for i in range(nSteps):
    # Perturb and check new structure:
    chains.randomRotation(1)
    bonds = chains.check()
    deltaBonds = oldBonds - bonds
    # Metropolis condition:
    if np.isnan(bonds):
        keep = False
    elif (deltaBonds <= 0) or (np.random.rand() < np.exp(-chains.beta*deltaBonds)):
        keep = True
    else:
        keep = False
    # Keep or restore:
    if keep:
        goodSteps += 1
        oldBonds = bonds
    else:
        chains.restorePrevious()
        bonds = oldBonds
        badSteps += 1
    # output
    i_ = i + 1
    if (i_ % outputFreq == 0) or (i_ == nSteps):
        t = time.time()
        print '   step %d/%d: %.1fs, %.1fs remaining'%(i_, nSteps, t-t0, float(nSteps-i_)*(t-t0)/i_)
        chains.dump(append=True)

# fix a nice vmd file for this simulation
fout = open(outputFile + '.vmd', 'w')
fin = open('base.vmd', 'r')
for line in fin:
    fout.write(line.replace('_FILENAME_', outputFile + '.pdb'))
if not surface:
    fout.write('\nsource align_traj.tcl \nalign_traj_on_itself 0 "all" 0\n')
fout.close()
fin.close()

print '\nSimulation done in %.1fs'%(t-t0) 
print 'Trajectory written to %s.pdb.'%outputFile
print 'Visualization state saved as %s.vmd, do "vmd -e %s.vmd" to look'%((outputFile,)*2)
print 'Acceptance rate: %.1f%%'%(100*float(goodSteps)/nSteps)
