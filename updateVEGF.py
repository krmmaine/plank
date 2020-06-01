from numpy import zeros
from math import cos
from math import pi

from xCoordinate import xCoordinate


def updateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength):
    k1 = 52923.07639
    k2 = 694.4444444
    k35 = 0.0129462015
    k21 = 10
    m0 = 12
    relax = 1.45
    k33 = 0.00036
    v = zeros((ySubstrate, xSteps))

    # capillary
    for x in range(xSteps-1):
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x+1]) / 2
        wallVEGFatJ = (vegf[1][x] + vegf[1][x+1]) / 2
        VEGFdiff = wallVEGFatJ - vegf[0][x]
        if VEGFdiff < 0:
            VEGFdiff = 0
        vegfOld[0][x] = vegf[0][x]
        vegf[0][x] = vegf[0][x] + k * (-k1 * vegf[0][x] * density / (1 + vegf[0][x]) + k2 * VEGFdiff)   # EQ 46
        if vegf[0][x] < 0:
            vegf[0][x] = 0

    # initialize v: guess
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                v[y][x] = vegf[y][x]
        else:
            for x in range(xSteps):
                v[y][x] = vegf[y][x]

    intol = 0
    while intol == 0:
        intol = 1

        # boundary at maximum y. includes source
        vOld = v[ySubstrate-1][0]
        v[ySubstrate-1][0] = v[ySubstrate-1][1]
        if v[ySubstrate-1][0]- vOld > tolerance or v[ySubstrate-1][0] - vOld < -tolerance:
            intol = 0
        for x in range(1, xSteps-2, 1):
            vOld = v[ySubstrate-1][x]
            # EQ 65 page 179
            v[ySubstrate-1][x] = k35 * h * ((1-cos(2 * pi * xCoordinate(x, ySubstrate-1, xSteps, xLength)))**m0) \
                                    + v[ySubstrate-3][x]
            if v[ySubstrate-1][x] - vOld > tolerance or v[ySubstrate-1][x] - vOld < -tolerance:
                intol = 0

        vOld = v[ySubstrate-1][xSteps-2]
        v[ySubstrate-1][xSteps-2] = v[ySubstrate-1][xSteps-3]
        if v[ySubstrate-1][xSteps-2] - vOld > tolerance or v[ySubstrate-1][xSteps-2] - vOld < -tolerance:
            intol = 0
            
        # boundary at maximum y - 1
        vOld = v[ySubstrate - 2][0]
        v[ySubstrate - 2][0] = v[ySubstrate - 2][1]
        if v[ySubstrate - 2][0] - vOld > tolerance or v[ySubstrate - 2][0] - vOld < -tolerance:
            intol = 0
        for x in range(1, xSteps-1, 1):
            vOld = v[ySubstrate - 2][x]
            v[ySubstrate - 2][x] = k35 * h * ((1-cos(2 * pi * xCoordinate(x, ySubstrate-2, xSteps, xLength))) ** m0)\
                                   + v[ySubstrate - 4][x]
            if v[ySubstrate - 2][x] - vOld > tolerance or v[ySubstrate - 2][x] - vOld < -tolerance:
                intol = 0
            vOld = v[ySubstrate - 2][xSteps - 1]
            v[ySubstrate - 2][xSteps - 1] = v[ySubstrate - 2][xSteps - 2]
            if v[ySubstrate - 2][xSteps - 1] - vOld > tolerance or v[ySubstrate - 2][xSteps - 1] - vOld < -tolerance:
                intol = 0

        # Cycle through interior rows
        for y in range(ySubstrate-3, 2, -1):
            # First meshpoint on row requires BC at x=0
            vOld = v[y][0]
            v[y][0] = v[y][1]
            if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                intol=0

            if y % 2 == 0:
                for x in range(1, xSteps-2, 1):
                    density = densityScale ** 2 * (occupiedOld[y//2][x]+occupiedOld[y//2][x+1])/2
                    vOld = v[y][x]
                    v[y][x] = relax / (h*h + 2*k21*k) * (0.5*k21*k * (v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x] +
                                                                  vegf[y][x+1]+vegf[y][x-1]+vegf[y+2][x]+vegf[y-2][x])
                                                    + (h*h - 2*k21*k - h*h*k*k1*density / (1+vegf[y][x])) * vegf[y][x])\
                            + (1-relax) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol = 0

                vOld = v[y][xSteps-2]
                v[y][xSteps-2] = v[y][xSteps-3]
                if v[y][xSteps-2] - vOld > tolerance or v[y][xSteps-2] - vOld < -tolerance:
                    intol=0
            else:
                for x in range(1, xSteps-1, 1):
                    density = densityScale**2 * (occupiedOld[(y-1)//2][x] + occupiedOld[(y+1)//2][x]) / 2
                    vOld = v[y][x]
                    v[y][x] = relax/(h*h + 2*k21*k) * (0.5*k21*k*(v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x]
                                                            + vegf[y][x+1] + vegf[y][x-1] + vegf[y+2][x] + vegf[y-2][x])
                                                       + (h*h - 2*k21*k - h*h*k*k1*density / (1+vegf[y][x]))*vegf[y][x])\
                              + (1-relax) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol=0
                vOld = v[y][xSteps-1]
                v[y][xSteps-1] = v[y][xSteps-2]
                if v[y][xSteps-1] - vOld > tolerance or v[y][xSteps-1] - vOld < -tolerance:
                    intol = 0

        # Boundary row at y = 2
        vOld = v[2][0]
        v[2][0] = v[2][1]
        if v[2][0] - vOld > tolerance or v[2][0] - vOld < -tolerance:
            intol = 0
        for x in range(1,xSteps-2,1):
            vOld = v[2][x]
            v[2][x] = 1 / (k33 + h) * (k33 * v[4][x] + h*vegfOld[0][x])
            if v[2][x] - vOld > tolerance or v[2][x] - vOld < -tolerance:
                intol = 0
        vOld = v[2][xSteps-2]
        v[2][xSteps-2] = v[2][xSteps-3]
        if v[2][xSteps-2] - vOld > tolerance or v[2][xSteps-2] - vOld < -tolerance:
            intol = 0

        # Boundary row at y = 1: capillary wall
        vOld = v[1][0]
        v[1][0] = v[1][1]
        if v[1][0] - vOld > tolerance or v[1][0] - vOld < -tolerance:
            intol = 0
        for x in range(1, xSteps-1, 1):
            vOld = v[1][x]
            v[1][x] = 1 / (k33 + h) * (k33 * v[3][x] + h*(vegfOld[0][x-1] + vegfOld[0][x])/2)
            if v[1][x] - vOld > tolerance or v[1][x] - vOld < -tolerance:
                intol = 0
        vOld = v[1][xSteps-1]
        v[1][xSteps-1] = v[1][xSteps-2]
        if v[1][xSteps-1] - vOld > tolerance or v[1][xSteps-1] - vOld < -tolerance:
            intol = 0

    # Cycle through VEGF meshpoints and set VEGF at time step j+1
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0
        else:
            for x in range(xSteps):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

    return
