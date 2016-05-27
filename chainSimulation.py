#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import Chains
import time
import sys

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
    print "   -ramps <n>:       ramp up beta from 0 <n> times for simulated annealing, " + b + "where <n> = 0 gives constant beta as specified by" + b + "-beta. After each ramp, beta is held constant for the" + b + "equivalent of a ramp duration. (default 0)"
    print "   -outputFile <s>:  the base name of the outputfile (default: 'out')"
    print "   -outputFreq <n>:  produces output every <n> simulation steps (default: 10)"
    print "   -steps <n>:       runs a simulation of <n> steps, including steps where the" + b + "change is rejected (default 1000)"
    print "   -append:          if present, continues the existing simulation specified by" + b + "-outputFile, which must match the other settings" + b + "(default no)"
    print ""
    exit()

# parse arguments
if '-help' in sys.argv: usage()
surface = '-surface' in sys.argv
append = '-append' in sys.argv
def parse(argv, key, default):
    if key in argv:
        return argv[argv.index(key) + 1]
    else:
        return default
number     = int(parse(sys.argv, '-number', 1))
length     = int(parse(sys.argv, '-length', 50))
box        = float(parse(sys.argv, '-box', 10))
maxAngle   = float(parse(sys.argv, '-maxAngle', 90)) * np.pi / 180.0
beta       = float(parse(sys.argv, '-beta', 0))
ramps      = int(parse(sys.argv, '-ramps', 0))
outputFile = parse(sys.argv, '-outputFile', 'out')
outputFreq = int(parse(sys.argv, '-outputFreq', 10))
nSteps     = int(parse(sys.argv, '-steps', 1000))
if (number > 1) and not surface:
    print "Simulating more than one chain without a surface makes no sense!\n"
    exit()

# start the timer and initialize the Chains instance
t0 = time.time()
initialConf = {False: None, True: outputFile+'.pdb'}[append]
chains = Chains.Chains(number=number, length=length, box=box, maxAngle=maxAngle, beta=beta, surface=surface, outFile=outputFile+'.pdb', initialConf=initialConf)

# set up the output file
if not append:
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
    # Find time-dependent beta for simulated annealing:
    if ramps == 0:
        beta_ = chains.beta
    else:
        beta_ = chains.beta
        if not (2 * i / (nSteps / (ramps))) % 2:
            beta_ = chains.beta * (2 * i % (nSteps / ramps)) / (nSteps / ramps)
    # Metropolis condition:
    if np.isnan(bonds):
        keep = False
    elif (deltaBonds <= 0) or (np.random.rand() < np.exp(-beta_*deltaBonds)):
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
        betaText  = ''
        if ramps:
            betaText = ', beta=%.2f'%beta_
        remaining = float(nSteps-i_)*((t-t0)%60)/i_
        timeString = {True:'%.1fs'%remaining, False:'%.0fmin'%(remaining//60,)}[remaining<300]
        print '   step %d/%d: %.1fs, %s remaining'%(i_, nSteps, t-t0, timeString) + betaText
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

t = time.time()
timeText = {True: '%.1fs'%(t-t0), False: '%.0fmin'%((t-t0)/60)}[t-t0 < 300]
print '\nSimulation done in ' + timeText
print 'Trajectory written to %s.pdb'%outputFile
print 'Visualization state saved as %s.vmd, do "vmd -e %s.vmd" to look'%((outputFile,)*2)
if nSteps > 0:
    print 'Acceptance rate: %.1f%%'%(100*float(goodSteps)/nSteps)
