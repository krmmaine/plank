from math import exp
from random import random
from heaviside import heaviside
from pStay import pStay
from pMove import pMove
from move import move
from updatePro import updatePro
from updateVEGF import updateVEGF
from newUpdateVEGF import newUpdateVEGF
from updateVEGFfib import updateVEGFfib
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
               movement, maxCell, birthTime, divideTime):
    # v0 = 0.04
    v0 = 0.2
    Dv = 3.6 * (10 ** -5)

    densityMax = 6
    k6 = 0.012048193
    k18 = 0.55555555
    k20 = 0.0496031736
    k25 = 5736.899771
    k26 = .00001859
    m1 = 2
    # child = 0.25      # 12 hours iteration 5400
    # child = 0.125     # 6 hours iteration 2700
    child = 0.0625      # 3 hours iteration 1350
    fibThreshold = 0.6
    workspace = zeros((ySteps, xSteps))
    file = open("tracking_backtracks.txt", "w")
    file2 = open("prolif_death.txt", "w")
    file3 = open("parameters.txt", "w")

    file3.write("PARAMETERS\n")
    file3.write("# hours until branching: " + str(48*child) + "\n")
    file3.write("VEGF flux term: " + str(v0) + "\n")
    file3.write("VEGF concentration: " + str(v0*Dv) + "\n")
    prolif = 0
    anastomosis = 0

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
            # not including dividing and death right now because it doesn't work with anastomosis
            #'''
            # DETERMINE IF CELL HAS DIVIDED OR DIED
            #'''
            if deathTime[cell] == numTimeSteps - 1 and y > 0:
                prolif = 1
                # add statement to test if dead
                # cell dies/leaves simulation if it reaches the tumour
                if y == ySteps - 1:
                    deathTime[cell] = time
                    occupied[y][x] -= 1

                # calculate average protease values at time and time-1
                # surrounding mesh points usually equals 4 unless it is a boundary condition
                proMinus0 = 0
                proMinus1 = 0
                count = 0

                #LEFT
                if x > 0:
                    proMinus0 = proMinus0 + pro[2*y][x-1] / (1 + k6 * fib[2*y][x-1])
                    proMinus1 = proMinus1 + proOld[2*y][x-1] / (1 + k6 * fibOld[2*y][x-1])
                    count += 1

                #RIGHT
                if x < xSteps - 1:
                    proMinus0 = proMinus0 + pro[2*y][x] / (1 + k6 * fib[2*y][x])
                    proMinus1 = proMinus1 + proOld[2*y][x] / (1 + k6 * fibOld[2*y][x])
                    count += 1

                #UP
                if y > 0:
                    proMinus0 = proMinus0 + pro[2*y-1][x] / (1 + k6 * fib[2*y-1][x])
                    proMinus1 = proMinus1 + proOld[2*y-1][x] / (1 + k6 * fibOld[2*y-1][x])
                    count += 1

                #DOWN
                if y < ySteps - 1:
                    proMinus0 = proMinus0 + pro[2*y+1][x-1] / (1 + k6 * fib[2*y+1][x-1])
                    proMinus1 = proMinus1 + proOld[2*y+1][x-1] / (1 + k6 * fibOld[2*y+1][x-1])
                    count += 1

                proMinus0 = proMinus0 / count
                proMinus1 = proMinus1 / count


                # logistic term is always 0 becauce denScale times occupied is always greater than densitymax
                if densityScale * occupied[y][x] > densityMax:
                    logistic = 0
                else:
                    logistic = 1 - densityScale * occupied[y][x] / densityMax

                G = k25 * (exp(-k26*proMinus0**m1)*(1-k26*m1*proMinus0**m1)) / (1+k25*proMinus0*exp(-k26*proMinus0**m1))

                proDependent = G * (proMinus0 - proMinus1) / k

                if proDependent >= 0:
                    # haven't figured out logistic term yet
                    # divideProb = (k * k18 + G * logistic * (proMinus0 - proMinus1)) * \
                    #             heaviside(time - divideTime[cell] - child * numTimeSteps)
                    divideProb = (k * k18 + G * (proMinus0 - proMinus1)) * \
                                 heaviside(time - divideTime[cell] - child * numTimeSteps)
                    deathProb = k * k20
                else:
                    divideProb = k * k18 * heaviside(time - divideTime[cell] - child * numTimeSteps)
                    deathProb = k * k20 - G * (proMinus0 - proMinus1)

                # not sure if I need this part
                #if divideProb < 0:
                #    divideProb = 0
                #if divideProb > 1 - deathProb:
                #    divideProb = 1

                randomNum = random()

                if randomNum < deathProb:
                    file2.write("\n\nTIME: " + str(time) + "\n")
                    file2.write("CELL: " + str(cell) + "\n")
                    file2.write("cell position: y = " + str(y) + " x = " + str(x) + "\n")
                    file2.write("average protease at time - 1 : " + str(proMinus1) + "\n")
                    file2.write("average protease at time: " + str(proMinus0) + "\n")
                    #file2.write("logistic growth term: " + str(logistic) + "\n")
                    file2.write("Protease dependent term: " + str(proDependent) + "\n")
                    file2.write("Probability of division: " + str(divideProb) + "\n")
                    file2.write("Probability of death: " + str(deathProb) + "\n")
                    file2.write("Probability of doing nothing: " + str(1 - deathProb - divideProb) + "\n")
                    file2.write("random number: " + str(randomNum) + "\n")
                    file2.write("G: " + str(G) + "\n")
                    file2.write("CELL DIED")
                    deathTime[cell] = time
                    occupied[y][x] -= 1
                    createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time)

                elif randomNum < deathProb + divideProb:
                    if totNumCells >= maxCell:
                        print("too many cells. DO NOT DIVIDE")
                    elif totNumCells < maxCell:
                        file2.write("\n\nTIME: " + str(time) + "\n")
                        file2.write("CELL: " + str(cell) + "\n")
                        file2.write("cell position: y = " + str(y) + " x = " + str(x) + "\n")
                        file2.write("average protease at time - 1 : " + str(proMinus1) + "\n")
                        file2.write("average protease at time: " + str(proMinus0) + "\n")
                        #file2.write("logistic growth term: " + str(logistic) + "\n")
                        file2.write("Protease dependent term: " + str(proDependent) + "\n")
                        file2.write("Probability of division: " + str(divideProb) + "\n")
                        file2.write("Probability of death: " + str(deathProb) + "\n")
                        file2.write("Probability of doing nothing: " + str(1 - deathProb - divideProb) + "\n")
                        file2.write("random number: " + str(randomNum) + "\n")
                        file2.write("G: " + str(G) + "\n")
                        file2.write("CELL DIVIDED")
                        divided = 1

                        if movement[len(movement)-1] == 0 or movement[len(movement)-1] == 1:
                            occupied[y][x + 1] += 1     # add cell to the right
                            occupied[y][x] -= 1         # subtract cell from current location
                            occupied[y][x - 1] += 1     # add cell to the left

                            createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time)

                            xPos[totNumCells][time+1] = x + 1           # new cell is to the right
                            yPos[totNumCells][time+1] = y

                            workspace[yPos[cell][time]][xPos[cell][time]] = 0

                            xPos[cell][time] = x - 1                # cell that divided is left
                            xPos[cell][time+1] = x - 1              # predicted location. if it stays it will be here.
                            x = xPos[cell][time]                    # reset x because position has changed

                        if movement[len(movement) - 1] == 2 or movement[len(movement) - 1] == 3:
                            occupied[y - 1][x] += 1  # add cell up
                            occupied[y][x] -= 1  # subtract cell from current location
                            occupied[y + 1][x] += 1  # add cell down

                            createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time)

                            xPos[totNumCells][time + 1] = x
                            yPos[totNumCells][time + 1] = y - 1     # new cell is up

                            workspace[yPos[cell][time]][xPos[cell][time]] = 0

                            yPos[cell][time] = y + 1                # cell that divided is down
                            yPos[cell][time + 1] = y + 1            # predicted location. if it stays it will be here.
                            y = yPos[cell][time]                    # reset y because position has changed

                        workspace[yPos[cell][time]][xPos[cell][time]] = (cell+1)*3
                        workspace[yPos[totNumCells][time + 1]][xPos[totNumCells][time + 1]] = (totNumCells+1)*3

                        deathTime[totNumCells] = numTimeSteps - 1
                        birthTime[totNumCells] = time
                        divideTime[cell] = time
                        divideTime[totNumCells] = time
                        totNumCells += 1
                        createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time+1)
            
            #'''

            # DETERMINE IF/WHERE THE CELL MOVES
            if deathTime[cell] == numTimeSteps - 1:
                stay = pStay(y, lamda, k)
                left, T = pMove(x, y, 0, pro, fib, vegf, xSteps, ySteps, lamda, k, movement)
                right, T = pMove(x, y, 1, pro, fib, vegf, xSteps, ySteps, lamda, k, movement)
                up, T = pMove(x, y, 2, pro, fib, vegf, xSteps, ySteps, lamda, k, movement)
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
                file = move(cell, time, stay, left, right, up, rand, yPos, xPos, occupied, fib, vegf, pro, movement, file, T, xSteps, ySteps)

                # ANASTOMOSIS
                '''
                if yPos[cell][time + 1] > 0:
                    anastomosis = 1
                    if workspace[yPos[cell][time + 1]][xPos[cell][time + 1]] != 0 and \
                            workspace[yPos[cell][time + 1]][xPos[cell][time + 1]] != (cell+1)*2 and \
                            workspace[yPos[cell][time + 1]][xPos[cell][time + 1]] != (cell+1)*3:
                        print("cell ran into another capillary")
                        deathTime[cell] = time + 1
                        occupied[yPos[cell][time + 1]][xPos[cell][time + 1]] -= 1

                workspace[yPos[cell][time]][xPos[cell][time]] = (cell+1)*2
                workspace[yPos[cell][time + 1]][xPos[cell][time + 1]] = (cell+1)*3
                '''


                # workspace without anastomosis
                workspace[yPos[cell][time]][xPos[cell][time]] = 1
                workspace[yPos[cell][time+1]][xPos[cell][time+1]] = 5


        total = 0
        for cell in range(totNumCells):
            if deathTime[cell] != numTimeSteps - 1:
                total += 1
        if total == totNumCells:
            print("ALL OF THE CELLS ARE DEAD")
            break

        updateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength, v0, Dv)
        # vegf = newUpdateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength, v0, Dv)
        updateFib(ySubstrate, xSteps, densityScale, occupiedOld, fib, fibOld, k, pro, tolerance, h)
        # updateVEGFfib(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, fib, fibOld, pro, proOld, k,
        #               tolerance, h, xLength, v0, Dv)
        updatePro(ySubstrate, xSteps, densityScale, occupiedOld, pro, proOld, k, vegfOld)

        for cell in range(totNumCells):
            if yPos[cell][time] > 0:
                fib[yPos[cell][time]][xPos[cell][time]] = 1

        print("time = " + str(time))


        if time % 500 == 0:
            createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time)


    if prolif == 1:
        file3.write("Branching and Death: YES\n")
    else:
        file3.write("Branching and Death: NO\n")

    if anastomosis == 1:
        file3.write("Anastomosis: YES\n")
    else:
        file3.write("Anastomosis: NO\n")

    return
