from numpy import zeros
from numpy import ones

from timePerStep import timePerStep
from init import init
from stepSize import stepSize
from lamda import lamda
from simulation import simulation

import time


def main():
    start = time.time()
    print(start)
    yLength = 0.5
    xLength = 1
    xSteps = 201
    ySteps = int(xSteps * (yLength/xLength) + 0.5)  # need to add the 0.5 so that it rounds to the correct number
    xSubstrate = xSteps * 2 - 1
    ySubstrate = ySteps * 2 - 1
    # totNumCells must be less than or equal to xStep so there will be a place to put all the cells
    totNumCells = 2
    densityScale = totNumCells/xSteps
    maxCell = 100
    numTimeSteps = 21600
    totalTime = 0.06912
    xPos = zeros((maxCell, numTimeSteps), dtype=int)
    yPos = zeros((maxCell, numTimeSteps), dtype=int)
    tolerance = 0.000001
    deathTime = zeros(maxCell, dtype=int)

    vegf = zeros((ySubstrate, xSteps))
    pro = zeros((ySubstrate, xSteps))
    fib = ones((ySubstrate, xSteps))
    vegfOld = zeros((ySubstrate, xSteps))
    proOld = zeros((ySubstrate, xSteps))
    fibOld = ones((ySubstrate, xSteps))

    occupied = zeros((ySteps, xSteps), dtype=int)
    occupiedOld = zeros((ySteps, xSteps), dtype=int)

    k = timePerStep(totalTime, numTimeSteps)
    h = stepSize(xLength, xSteps)
    lam = lamda(xLength, h)

    init(xSteps, totNumCells, xPos, yPos, occupied, deathTime, numTimeSteps)
    simulation(numTimeSteps, xSteps, ySteps, occupied, occupiedOld, totNumCells, xPos, yPos, deathTime, pro, proOld,
               densityScale, lam, k, fib, vegf, ySubstrate, vegfOld, tolerance, h, xLength, fibOld)
    print("finished yay!!")
    print(time.time())
    print(time.time()-start)

    return


'''
import cProfile
import re
cProfile.run('re.compile("foo|bar")', 'restats')

import pstats
p = pstats.Stats('restats')
p.strip_dirs().sort_stats(-1).print_stats()

p.sort_stats('name')
p.print_stats()

print("profiling done")
'''

main()
