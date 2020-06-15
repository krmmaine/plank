from math import exp
from random import random
from heaviside import heaviside
from pStay import pStay
from pMove import pMove
from move import move
from updatePro import updatePro
from updateVEGF import updateVEGF
from updateFib import updateFib
from numpy import zeros

from scipy import spatial
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def simulation(numTimeSteps, xSteps, ySteps, occupied, occupiedOld, totNumCells, xPos, yPos, deathTime, pro, proOld, densityScale, lamda, k, fib, vegf, ySubstrate, vegfOld, tolerance, h, xLength, fibOld):
    densityMax = 6
    k25 = 5736.899771
    k26 = .00001859
    m1 = 2
    fibThreshold = 0.6
    workspace = zeros((ySteps, xSteps))

    # Cycle through time steps
    for time in range(numTimeSteps-1):
        # Copy occupied into occupiedOld
        for x in range(xSteps):
            for y in range(ySteps):
                occupiedOld[y][x] = occupied[y][x]
        # Cycle through cells
        for cell in range(totNumCells):
            # At first assume no movement
            x = xPos[cell][time]
            y = yPos[cell][time]
            xPos[cell][time+1] = x
            yPos[cell][time+1] = y

            # If cell has left the capillary
            # DETERMINE IF CELL HAS DIVIDED OR DIED
            if deathTime[cell] == numTimeSteps - 1 and y > 0:
                # cell dies/leaves simulation if it reaches the tumour
                if y == ySteps - 1:
                    deathTime[cell] = time
                    occupied[y][x] -= 1
                # Calculate average protease at surrounding substrate meshpoints at time and time minus 1
                proAtTime = 0
                proAtTimeMinus1 = 0
                count = 0
                # Protease to the left
                if x > 0:
                    proAtTime = proAtTime + pro[2*y][x-1]
                    proAtTimeMinus1 = proAtTimeMinus1 + proOld[2*y][x-1]
                    count += 1
                # Protease to the right
                if x < xSteps - 1:
                    proAtTime = proAtTime + pro[2*y][x]
                    proAtTimeMinus1 = proAtTimeMinus1 + proOld[2*y][x]
                    count += 1
                # Protease below
                if y > 0:
                    proAtTime = proAtTime + pro[2*y-1][x]
                    proAtTimeMinus1 = proAtTimeMinus1 + proOld[2*y-1][x]
                    count += 1
                # Protease above
                if y < ySteps - 1:
                    proAtTime = proAtTime + pro[2*y+1][x]
                    proAtTimeMinus1 = proAtTimeMinus1 + proOld[2*y+1][x]
                    count += 1

                proAtTime = proAtTime/count
                proAtTimeMinus1 = proAtTimeMinus1/count

                # If protease level has gone down then there is no contribution to proliferation
                if proAtTimeMinus1 > proAtTime:
                    proAtTimeMinus1 = proAtTime
                if densityScale * occupied[y][xSteps] > densityMax:
                    logistic = 0
                else:
                    logistic = 1 - densityScale * occupied[y][x] / densityMax

                G = k25 * exp(-k26 * proAtTime ** m1) * (1 - k26 * m1 * proAtTime ** m1) / \
                    (1 + k25 * proAtTime * exp(-k26 * proAtTime ** m1))
                if G < 0:
                    G = 0

                # NEED TO FIGURE OUT PROLIFERATION AND CELL DEATH

                # DETERMINE IF/WHERE THE CELL MOVES

                if deathTime[cell] == numTimeSteps - 1:
                    stay = pStay(y, lamda, k)
                    left = pMove(x, y, 0, pro, fib, vegf, xSteps, ySteps, lamda, k)
                    right = pMove(x, y, 1, pro, fib, vegf, xSteps, ySteps, lamda, k)
                    down = pMove(x, y, 2, pro, fib, vegf, xSteps, ySteps, lamda, k)
                    rand = random()
                    # Check if cell can escape the capillary
                    if y == 0:
                        if x == 0:
                            fibcap = fib[0][0]
                        elif x == xSteps - 1:
                            fibcap = fib[0][xSteps-2]
                        else:
                            fibcap = (fib[0][x-1] + fib[0][x]) / 2
                        if fibcap < fibThreshold:
                            rand = 2
                    move(cell, time, stay, left, right, down, rand, yPos, xPos, occupied)
                    workspace[xPos[cell][time]][yPos[cell][time]] = 1000

        updateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength)
        updateFib(ySubstrate, xSteps, densityScale, occupiedOld, fib, fibOld, k, pro, tolerance, h)
        updatePro(ySubstrate, xSteps, densityScale, occupiedOld, pro, proOld, k, vegfOld)

        if time % 1000 == 0:
            plt.imshow(workspace)
            cm.get_cmap("jet")
            plt.show()

            fig = plt.figure()
            ax = plt.gca(projection='3d')

            X = np.arange(0, xSteps, 1)
            Y = np.arange(0, ySubstrate, 1)
            X, Y = np.meshgrid(X, Y)

            Z = vegf

            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                   cmap='viridis', edgecolor='none')
            ax.set_title('vegf')
            plt.show()

            fig = plt.figure()
            ax = plt.gca(projection='3d')

            X = np.arange(0, xSteps, 1)
            Y = np.arange(0, ySubstrate, 1)
            X, Y = np.meshgrid(X, Y)

            Z = fib

            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                   cmap='viridis', edgecolor='none')
            ax.set_title('fibronectin')
            plt.show()

            fig = plt.figure()
            ax = plt.gca(projection='3d')

            X = np.arange(0, xSteps, 1)
            Y = np.arange(0, ySubstrate, 1)
            X, Y = np.meshgrid(X, Y)

            Z = pro

            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                   cmap='viridis', edgecolor='none')
            ax.set_title('protease')
            plt.show()
    print("simulation loop complete. check variables")

    return
