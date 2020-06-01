from numpy import zeros
from numpy import ones

from timePerStep import timePerStep
from init import init
from stepSize import stepSize
from lamda import lamda
from simulation import simulation


def main():
    yLength = 0.5
    xLength = 1
    xSteps = 200
    ySteps = int(xSteps * (yLength/xLength))
    xSubstrate = xSteps * 2 - 1
    ySubstrate = ySteps * 2 - 1
    # totNumCells must be less than or equal to xStep so there will be a place to put all the cells
    totNumCells = 5
    densityScale = totNumCells/xSteps
    maxCell = 100
    numTimeSteps = 21600
    totalTime = 0.06912
    xPos = zeros((maxCell, numTimeSteps))
    yPos = zeros((maxCell, numTimeSteps))
    tolerance = 0.000001
    deathTime = zeros(maxCell)

    vegf = zeros((ySubstrate, xSteps))
    pro = zeros((ySubstrate, xSteps))
    fib = ones((ySubstrate, xSteps))
    vegfOld = zeros((ySubstrate, xSteps))
    proOld = zeros((ySubstrate, xSteps))
    fibOld = ones((ySubstrate, xSteps))

    occupied = zeros((ySteps, xSteps))
    occupiedOld = zeros((ySteps, xSteps))

    k = timePerStep(totalTime, numTimeSteps)
    h = stepSize(xLength, xSteps)
    lam = lamda(xLength, h)

    init(xSteps, totNumCells, xPos, yPos, occupied, deathTime, numTimeSteps)
    simulation(numTimeSteps, xSteps, ySteps, occupied, occupiedOld, totNumCells, xPos, yPos, deathTime, pro, proOld,
               densityScale, lam, k, fib, vegf, ySubstrate, vegfOld, tolerance, h, xLength, fibOld)
    print("finished yay!!")

    return


main()
