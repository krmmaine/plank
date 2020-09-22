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
from graph import createGraph

from scipy import spatial
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import matplotlib.backends.backend_pdf


def simulation(numTimeSteps, xSteps, ySteps, occupied, occupiedOld, totNumCells, xPos, yPos, deathTime, pro, proOld,
               densityScale, lamda, k, fib, vegf, ySubstrate, vegfOld, tolerance, h, xLength, fibOld, xVector, yVector,
               movement):
    densityMax = 6
    k25 = 5736.899771
    k26 = .00001859
    m1 = 2
    fibThreshold = 0.6
    workspace = zeros((ySteps, xSteps))
    file = open("tracking_backtracks.txt", "w")

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

            # DETERMINE IF/WHERE THE CELL MOVES

            if deathTime[cell] == numTimeSteps - 1:
                stay = pStay(y, lamda, k)
                left = pMove(x, y, 0, pro, fib, vegf, xSteps, ySteps, lamda, k)
                right = pMove(x, y, 1, pro, fib, vegf, xSteps, ySteps, lamda, k)
                up = pMove(x, y, 2, pro, fib, vegf, xSteps, ySteps, lamda, k)
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
                file = move(cell, time, stay, left, right, up, rand, yPos, xPos, occupied, fib, vegf, pro, movement, file)
                workspace[yPos[cell][time]][xPos[cell][time]] = 2
                workspace[yPos[cell][time+1]][xPos[cell][time+1]] = 5

        updateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength)
        updateFib(ySubstrate, xSteps, densityScale, occupiedOld, fib, fibOld, k, pro, tolerance, h)
        updatePro(ySubstrate, xSteps, densityScale, occupiedOld, pro, proOld, k, vegfOld)

        print("time = " + str(time))

        if time % 2000 == 0:
            createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time)

    return
