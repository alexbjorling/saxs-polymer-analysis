import numpy as np
#import multiprocessing
#import time
import copy
import sys

class Chains(object):
    def __init__(self, number=0, length=0, box=0, maxAngle=np.pi/2, beta=0.0, surface=False, outFile='out.pdb', initialConf=None):
        # initiate a set of <number> chains of length <length>, spread randomly over a square surface of side <box>
        if number > 0:
            self.number = number
            self.length = length
            self.box = box
            self.maxAngle = maxAngle
            self.beta = beta
            self.surface = surface
            self.outFile = outFile
            bonds = np.nan
            # read the initial conformation from pdb
            if initialConf:
                print "Reading initial conformation from the last frame of " + initialConf + "."
                # go to the first line of the last frame
                fp = open(initialConf, 'r')
                i = 0
                for line in fp:
                    i += 1
                    if line.startswith('MODEL'):
                        linesBefore = i
                fp.seek(0)
                for i in range(linesBefore):
                    fp.readline()
                self.coords = []
                for i in range(number):
                    self.coords.append(np.ones((length, 3)))
                    for j in range(length):
                        line = fp.readline()
                        self.coords[-1][j,:] = np.array(map(float, [line[30:38], line[38:46], line[46:54]]))
                    fp.readline() # TER
                fp.close()
            # make a random one
            else:
                while np.isnan(bonds):
                    sys.stdout.write("Generating initial arrangement of chains... ")
                    self.coords = []
                    for i in range(number):
                        self.coords.append(np.ones((length, 3)))
                        self.coords[-1][:,0] = np.random.rand() * box
                        self.coords[-1][:,1] = np.random.rand() * box
                        self.coords[-1][:,2] = np.arange(.5, length)
                    bonds = self.check()
                    if np.isnan(bonds):
                        sys.stdout.write('failed. Trying again.\n')
                    else:
                        sys.stdout.write('success!\n')
                
    def randomRotation(self, number=1):
        for cycle in range(number):
            m = np.random.randint(self.number)
            n = np.random.randint(self.length-1)
            thetaMax = self.maxAngle # set the random angle to the same value as the maximum allowed angle - should be the right sort of size
            tX = (1-2*np.random.rand()) * thetaMax * (1 + self.surface)
            tY = (1-2*np.random.rand()) * thetaMax * (1 + self.surface)
            tZ = (1-2*np.random.rand()) * thetaMax * (1 + self.surface)
            rotX = np.array([[1, 0, 0], [0, np.cos(tX), -np.sin(tX)], [0, np.sin(tX), np.cos(tX)]])
            rotY = np.array([[np.cos(tY), 0, np.sin(tY)], [0, 1, 0], [-np.sin(tY), 0, np.cos(tY)]])
            rotZ = np.array([[np.cos(tZ), -np.sin(tZ), 0], [np.sin(tZ), np.cos(tZ), 0], [0, 0, 1]])
            rot = rotX.dot(rotY).dot(rotZ)
            origin = self.coords[m][n+1, :]
            self.oldCoords = copy.deepcopy(self.coords)
            for i in range(n+2, self.length):
                self.coords[m][i,:] = origin + np.dot(rot, self.coords[m][i,:]-origin)
                
    def restorePrevious(self):
        if self.oldCoords == None:
            raise RuntimeError('No saved state to restore')
        else:
            self.coords = self.oldCoords
            self.oldCoords = None
            
    def check(self):
        bonds = 0
        # got to loop over chains first, beads i and i+1 should interact fully if they are on different chains...
        for m in range(self.number):
            # check for angle violation on the current chain:
            for i in range(1, self.length-1):
                a = self.coords[m][i+1,:] - self.coords[m][i,:]
                b = self.coords[m][i,:] - self.coords[m][i-1,:]
                angle = np.arccos( min(1, np.dot(a,b) / (np.linalg.norm(a)*np.linalg.norm(b)) ))
                if angle > self.maxAngle:
                    return np.nan
            # check for overlap and surface violation on the current chain:
            for i in range(self.length):
                ri = self.coords[m][i]
                if self.surface and (i>0) and (ri[2] < 0.5):
                    return np.nan
                for j in range(i-1):
                    rj = self.coords[m][j]
                    if np.linalg.norm(ri-rj) < 1.0:
                        return np.nan
                    # count bonds
                    if np.linalg.norm(ri-rj) < 1.2:
                        bonds += 1
            # mn cross terms:
            for n in range(m):
                # check for overlap and bonds
                for i in range(self.length):
                    ri = self.coords[m][i]
                    for j in range(self.length):
                        rj = self.coords[n][j]
                        if np.linalg.norm(ri-rj) < 1.0:
                            return np.nan
                        # count bonds
                        if np.linalg.norm(ri-rj) < 1.2:
                            bonds += 1
        return bonds
            
    def dump(self, append=False):
        fp = open(self.outFile, {True:'a', False:'w'}[append])
        fp.write('MODEL \n')
        for j in range(self.number):
            for i in range(self.length):
                if i == 0: 
                    name = 'A'
                elif i == 1:
                    name = 'B'
                else:
                    name = 'C'
                i_ = j * self.length + i
                fp.write("ATOM    %3d    %s AAA A %3d    %8.3f%8.3f%8.3f\n"%(i_, name, j, self.coords[j][i,0], self.coords[j][i,1], self.coords[j][i,2]))
            fp.write('TER   \n')
        fp.write('ENDMDL\n')
        fp.close()
        
    def copy(self):
        return copy.deepcopy(self)