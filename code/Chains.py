"""
Contains the man class to represent chains, perturb their structures,
check for violations and count bonds.
"""

import numpy as np
import copy
import sys


class Chains(object):
    def __init__(self, number=0, length=0, box=0, maxAngle=np.pi / 2,
                 beta=0.0, surface=False, outFile='out.pdb', initialConf=None,
                 grid=None):
        """
        Make and initialize a set of <number> chains of length <length>,
        spread over a square surface of side <box>, MxN if <grid> is
        specified, otherwise randomly.
        """
        if number > 0:
            self.number = number
            self.length = length
            self.box = box
            self.maxAngle = maxAngle
            self.beta = beta
            self.surface = surface
            self.outFile = outFile
            bonds = None
            # read the initial conformation from pdb
            if initialConf:
                print("Reading initial conformation from the last frame of %s."
                      % initialConf)
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
                        self.coords[-1][j, :] = np.array(
                            map(float,
                                [line[30:38], line[38:46], line[46:54]])
                        )
                    fp.readline()  # TER
                fp.close()
            # make a random one
            else:
                while bonds is None:
                    sys.stdout.write(
                        "Generating initial arrangement of chains... ")
                    self.coords = []
                    if grid:
                        assert grid[0] * grid[1] >= number
                        dx = box / (grid[0] + 1)
                        x = np.linspace(0, 1, grid[0] + 1)[:-1] * box + dx / 2
                        dy = box / (grid[1] + 1)
                        y = np.linspace(0, 1, grid[1] + 1)[:-1] * box + dy / 2
                        X, Y = np.meshgrid(x, y)
                        X = X.flatten()
                        Y = Y.flatten()
                    else:
                        X = np.random.rand(number) * box
                        Y = np.random.rand(number) * box
                    for i in range(number):
                        self.coords.append(np.ones((length, 3)))
                        self.coords[-1][:, 0] = X[i]
                        self.coords[-1][:, 1] = Y[i]
                        self.coords[-1][:, 2] = np.arange(.5, length)
                    bonds, angles = self.check()
                    if bonds is None:
                        sys.stdout.write('failed. Trying again.\n')
                    else:
                        sys.stdout.write('success!\n')

    def randomRotation(self, number=1, size=1):
        for cycle in range(number):
            m = np.random.randint(self.number)
            n = np.random.randint(self.length - 1)
            # set the random angle to the same value as the maximum
            # allowed angle - should be the right sort of size:
            thetaMax = self.maxAngle * size
            tX = (1 - 2 * np.random.rand()) * thetaMax * (1 + self.surface)
            tY = (1 - 2 * np.random.rand()) * thetaMax * (1 + self.surface)
            tZ = (1 - 2 * np.random.rand()) * thetaMax * (1 + self.surface)
            rotX = np.array(
                [[1, 0, 0], [0, np.cos(tX), -np.sin(tX)],
                 [0, np.sin(tX), np.cos(tX)]]
            )
            rotY = np.array(
                [[np.cos(tY), 0, np.sin(tY)], [0, 1, 0],
                 [-np.sin(tY), 0, np.cos(tY)]]
            )
            rotZ = np.array(
                [[np.cos(tZ), -np.sin(tZ), 0],
                 [np.sin(tZ), np.cos(tZ), 0], [0, 0, 1]]
            )
            rot = rotX.dot(rotY).dot(rotZ)
            origin = self.coords[m][n + 1, :]
            self.oldCoords = copy.deepcopy(self.coords)
            for i in range(n + 2, self.length):
                self.coords[m][i, :] = (
                    origin + np.dot(rot, self.coords[m][i, :] - origin))

    def restorePrevious(self):
        if self.oldCoords is None:
            raise RuntimeError('No saved state to restore')
        else:
            self.coords = self.oldCoords
            self.oldCoords = None

    def check(self):
        bonds = 0
        angles = []
        # got to loop over chains first, beads i and i+1 should interact
        # fully if they are on different chains...
        for m in range(self.number):
            angles_ = []
            # check for angle violation on the current chain:
            for i in range(1, self.length - 1):
                a = self.coords[m][i + 1, :] - self.coords[m][i, :]
                b = self.coords[m][i, :] - self.coords[m][i - 1, :]
                angle = np.arccos(min(1, np.dot(a, b) / (np.linalg.norm(a)
                                  * np.linalg.norm(b))))
                if angle > self.maxAngle:
                    return (None, None)
                angles_.append(angle)
            # check for overlap and surface/box violation on the current chain:
            for i in range(self.length):
                ri = self.coords[m][i]
                if self.surface and (i > 0) and (ri[2] < 0.5):
                    return (None, None)
                if self.box:
                    if ((ri[0] < 0) or (ri[0] > self.box) or (ri[1] < 0) or (ri[1] > self.box)):
                        return (None, None)
                for j in range(i - 1):
                    rj = self.coords[m][j]
                    if np.linalg.norm(ri - rj) < 1.0:
                        return (None, None)
                    # count bonds
                    if np.linalg.norm(ri - rj) < 1.2:
                        bonds += 1
            # mn cross terms:
            for n in range(m):
                # check for overlap and bonds
                for i in range(self.length):
                    ri = self.coords[m][i]
                    for j in range(self.length):
                        rj = self.coords[n][j]
                        if np.linalg.norm(ri - rj) < 1.0:
                            return (None, None)
                        # count bonds
                        if np.linalg.norm(ri - rj) < 1.2:
                            bonds += 1
            angles.append(angles_)
        return bonds, np.array(angles)

    def debye(self, q):
        res = np.zeros_like(q)
        coords = np.array(self.coords).reshape((-1, 3))
        for i in range(coords.shape[0]):
            res[:] = res[:] + 1
            for j in range(i):
                r = np.linalg.norm(coords[i] - coords[j])
                res[0] += 2
                res[1:] = res[1:] + 2 * np.sin(q[1:] * r) / (q[1:] * r)
        return res / res[0]

    def dump(self, append=False):
        fp = open(self.outFile, {True: 'a', False: 'w'}[append])
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
                fp.write(
                    "ATOM    %3d    %s AAA A %3d    %8.3f%8.3f%8.3f\n"
                    % (i_, name, j,
                       self.coords[j][i, 0],
                       self.coords[j][i, 1],
                       self.coords[j][i, 2])
                )
            fp.write('TER   \n')
        fp.write('ENDMDL\n')
        fp.close()

    def copy(self):
        return copy.deepcopy(self)
