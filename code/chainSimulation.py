"""
Main script which runs simulations with command-line argument parameters.
"""

import numpy as np
import os
import Chains
import time
import sys


def usage():
    msg = (
        '\nUsage:\n\n'
        'chainSimulation.py\n\n'
        '  -surface:         if present, simulates the chains grafted on a\n'
        '                    hard but inert surface (default no)\n'
        '  -number <n>:      simulates a collection of n chains (default 1)\n'
        '  -length <n>:      makes coils of n beads, with contour lengths\n'
        '                    n - 1 (default 50)\n'
        '  -box <d>:         grafts the chains on a square of side length d\n'
        '                    (default 10)\n'
        '  -maxAngle <a>:    the maximum bend angle at each bead in degrees\n'
        '                    where a = 0.0 is a straight rod (default 90.0)\n'
        '  -beta <b>:        the effective bond strength b = energy / (kT),\n'
        '                    where the well width is 0.2 (default 0.0)\n'
        '  -ramps <n>:       ramp up beta from 0 n times for simulated\n'
        '                    annealing, where n = 0 gives constant beta as\n'
        '                    specified by -beta. After each ramp, beta is\n'
        '                    held constant for the equivalent of a ramp\n'
        '                    duration. (default 0)\n'
        '  -outputFile <s>:  the base name of the outputfile (default "out")\n'
        '  -outputFreq <n>:  produces output every n simulation steps \n'
        '                    (default: 10)\n'
        '  -steps <n>:       runs a simulation of n steps, including steps\n'
        '                    where the change is rejected (default 1000)\n'
        '  -append:          if present, continues the existing simulation\n'
        '                    specified by -outputFile, which must match the\n'
        '                    other settings (default no)\n'
    )
    print(msg)
    exit()


def parse(argv, key, default):
    if key in argv:
        return argv[argv.index(key) + 1]
    else:
        return default


# parse arguments
if '-help' in sys.argv:
    usage()
surface = '-surface' in sys.argv
append = '-append' in sys.argv
number = int(parse(sys.argv, '-number', 1))
length = int(parse(sys.argv, '-length', 50))
box = float(parse(sys.argv, '-box', 10))
maxAngle = float(parse(sys.argv, '-maxAngle', 90)) * np.pi / 180.0
beta = float(parse(sys.argv, '-beta', 0))
ramps = int(parse(sys.argv, '-ramps', 0))
outputFile = parse(sys.argv, '-outputFile', 'out')
outputFreq = int(parse(sys.argv, '-outputFreq', 10))
nSteps = int(parse(sys.argv, '-steps', 1000))
if (number > 1) and not surface:
    print("Simulating more than one chain without a surface makes no sense!\n")
    exit()

# start the timer and initialize the Chains instance
t0 = time.time()
initialConf = {False: None, True: outputFile + '.pdb'}[append]
chains = Chains.Chains(
    number=number, length=length, box=box, maxAngle=maxAngle, beta=beta,
    surface=surface, outFile=outputFile + '.pdb', initialConf=initialConf,
)

# set up the output
fp = open(outputFile + '.commandLine', 'w')
for s in sys.argv:
    fp.write(s + ' ')
fp.write('\n')
fp.close()
if not append:
    try:
        os.remove(outputFile + '.pdb')
    except:
        pass
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
    elif (deltaBonds <= 0) or (np.random.rand() < np.exp(-beta_ * deltaBonds)):
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
        betaText = ''
        if ramps:
            betaText = ', beta=%.2f' % beta_
        remaining = float(nSteps - i_) * ((t - t0)) / i_
        if remaining < 300:
            timeString = '%.1fs' % remaining
        else:
            timeString = '%.0fmin' % (remaining // 60,)
        print('   step %d/%d: %.1fs, %s remaining'
              % (i_, nSteps, t - t0, timeString) + betaText)
        chains.dump(append=True)

# fix a nice vmd file for this simulation
fout = open(outputFile + '.vmd', 'w')

fin = open(os.path.dirname(os.path.realpath(__file__)) + '/base.vmd', 'r')
for line in fin:
    fout.write(line.replace('_FILENAME_', outputFile + '.pdb'))
if not surface:
    fout.write('\nsource align_traj.tcl \nalign_traj_on_itself 0 "all" 0\n')
fout.close()
fin.close()

t = time.time()
if t - t0 < 300:
    timeText = '%.1fs' % (t - t0)
else:
    timeText = '%.0fmin' % ((t - t0) / 60)
print('\nSimulation done in ' + timeText)
print('Trajectory written to %s.pdb' % outputFile)
print('Visualization state saved as %s.vmd, do "vmd -e %s.vmd" to look'
      % ((outputFile,) * 2))
if nSteps > 0:
    print('Acceptance rate: %.1f%%' % (100 * float(goodSteps) / nSteps))
