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

# update VEGF concentration in capillary
    for x in range(xSteps - 1):
        # scaled by densityScale to get density instead of the number of cells
        # average density at cell meshpoints to the right and left of the substrate meshpoint
        # use occupiedOld because occupied has already been updated this timestep
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x + 1]) / 2
        # average VEGF concentration on capillary wall at time j
        wallVEGFatJ = (vegf[1][x] + vegf[1][x + 1]) / 2
        VEGFdiff = wallVEGFatJ - vegf[0][x]
        if VEGFdiff < 0:
            VEGFdiff = 0
        vegfOld[0][x] = vegf[0][x]
        # update VEGF concentration in capillary using equation 46
        vegf[0][x] = vegf[0][x] + k * (-k1 * vegf[0][x] * density / (1 + vegf[0][x]) + k2 * VEGFdiff)
        if vegf[0][x] < 0:
            vegf[0][x] = 0

# initialize v: value of the previous time step
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                v[y][x] = vegf[y][x]
        else:
            for x in range(xSteps):
                v[y][x] = vegf[y][x]

# Keep iterating until the value of v at each meshpoint changes by less than a certain tolerance
    intol = 0
    while intol == 0:
        intol = 1

    # update VEGF concentration at boundary at maximum y. includes source
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[ySubstrate - 1][0]
        v[ySubstrate - 1][0] = v[ySubstrate - 1][1]
        # if the change in any v and vOld is greater than the tolerance then this loop will continue to run
        if v[ySubstrate - 1][0] - vOld > tolerance or v[ySubstrate - 1][0] - vOld < -tolerance:
            intol = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[ySubstrate - 1][x]
            # use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary
            v[ySubstrate - 1][x] = k35 * h * ((1 - cos(2 * pi * xCoordinate(x, ySubstrate - 1, xSteps, xLength))) ** m0) \
                                   + v[ySubstrate - 3][x]
            if v[ySubstrate - 1][x] - vOld > tolerance or v[ySubstrate - 1][x] - vOld < -tolerance:
                intol = 0

        # using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[ySubstrate - 1][xSteps - 2]
        v[ySubstrate - 1][xSteps - 2] = v[ySubstrate - 1][
            xSteps - 3]  # minus 3 because row is even
        if v[ySubstrate - 1][xSteps - 2] - vOld > tolerance or v[ySubstrate - 1][xSteps - 2] - vOld < -tolerance:
            intol = 0

    # update VEGF concentration at boundary at maximum y - 1
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[ySubstrate - 2][0]
        v[ySubstrate - 2][0] = v[ySubstrate - 2][1]
        if v[ySubstrate - 2][0] - vOld > tolerance or v[ySubstrate - 2][0] - vOld < -tolerance:
            intol = 0

        for x in range(1, xSteps - 1, 1):
            vOld = v[ySubstrate - 2][x]
            # use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary - 1
            v[ySubstrate - 2][x] = k35 * h * ((1 - cos(2 * pi * xCoordinate(x, ySubstrate - 2, xSteps, xLength))) ** m0) \
                                   + v[ySubstrate - 4][x]  # -4 because substrate points are at half meshpoints
            if v[ySubstrate - 2][x] - vOld > tolerance or v[ySubstrate - 2][x] - vOld < -tolerance:
                intol = 0

        # using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[ySubstrate - 2][xSteps - 1]
        v[ySubstrate - 2][xSteps - 1] = v[ySubstrate - 2][xSteps - 2]
        if v[ySubstrate - 2][xSteps - 1] - vOld > tolerance or v[ySubstrate - 2][xSteps - 1] - vOld < -tolerance:
            intol = 0

    # Cycle through interior rows and calculate new VEGF concentration
        for y in range(ySubstrate - 3, 2, -1):
            # using equation 70 derivation on page 179 calculate vegf at x = 0
            vOld = v[y][0]
            v[y][0] = v[y][1]
            if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                intol = 0

            # if row is even number of substrate meshpoints in x is xsteps-1
            if y % 2 == 0:
                for x in range(1, xSteps - 2, 1):  # minus 2 because last meshpoint has a boundary condition
                    # densityScale is squared because density in ECM is per unit area not length
                    # average density of cell meshpoint to the right and left
                    density = densityScale ** 2 * (occupiedOld[y // 2][x] + occupiedOld[y // 2][x + 1]) / 2
                    vOld = v[y][x]
                    # Approximate equation 53 using the crank-nicolson method see derivation on page 178
                    # relax is the successive over-relaxation term
                    v[y][x] = relax / (h * h + 2 * k21 * k) * (
                                0.5 * k21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x] +
                                                 vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                + (h * h - 2 * k21 * k - h * h * k * k1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - relax) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol = 0

                # using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][xSteps - 2]
                v[y][xSteps - 2] = v[y][xSteps - 3]
                if v[y][xSteps - 2] - vOld > tolerance or v[y][xSteps - 2] - vOld < -tolerance:
                    intol = 0
            else:
                # if row is odd number of substrate meshpoints in x is nn
                for x in range(1, xSteps - 1, 1):
                    # average density of cell meshpoints above and below
                    # y // 2 because there are twice as many points in the y direction because substrate meshpoints are at 1/2
                    density = densityScale ** 2 * (occupiedOld[(y - 1) // 2][x] + occupiedOld[(y + 1) // 2][x]) / 2
                    vOld = v[y][x]
                    v[y][x] = relax / (h * h + 2 * k21 * k) * (
                                0.5 * k21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x]
                                                 + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                + (h * h - 2 * k21 * k - h * h * k * k1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - relax) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol = 0

                # using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][xSteps - 1]
                v[y][xSteps - 1] = v[y][xSteps - 2]
                if v[y][xSteps - 1] - vOld > tolerance or v[y][xSteps - 1] - vOld < -tolerance:
                    intol = 0

    # calculate VEGF concentration at boundary row y = 2
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[2][0]
        v[2][0] = v[2][1]
        if v[2][0] - vOld > tolerance or v[2][0] - vOld < -tolerance:
            intol = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[2][x]
            # use equation 61 see derivation on page 179 of paper and pg 137 of notes
            # use vegfOld because capillary values have already been updated
            v[2][x] = 1 / (k33 + h) * (k33 * v[4][x] + h * vegfOld[0][x])
            if v[2][x] - vOld > tolerance or v[2][x] - vOld < -tolerance:
                intol = 0

        # using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[2][xSteps - 2]
        v[2][xSteps - 2] = v[2][xSteps - 3]
        if v[2][xSteps - 2] - vOld > tolerance or v[2][xSteps - 2] - vOld < -tolerance:
            intol = 0

    # calculate VEGF concentratin at boundary row y = 1: capillary wall
        vOld = v[1][0]
        v[1][0] = v[1][1]
        if v[1][0] - vOld > tolerance or v[1][0] - vOld < -tolerance:
            intol = 0
        for x in range(1, xSteps - 1, 1):
            vOld = v[1][x]
            # Take average because there is no meshpoint in the capillary directly below
            v[1][x] = 1 / (k33 + h) * (k33 * v[3][x] + h * (vegfOld[0][x - 1] + vegfOld[0][x]) / 2)
            if v[1][x] - vOld > tolerance or v[1][x] - vOld < -tolerance:
                intol = 0
        # using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[1][xSteps - 1]
        v[1][xSteps - 1] = v[1][xSteps - 2]
        if v[1][xSteps - 1] - vOld > tolerance or v[1][xSteps - 1] - vOld < -tolerance:
            intol = 0

    # Cycle through VEGF meshpoints and set VEGF at time step j+1
    for y in range(1, ySubstrate, 1):
        # if y is even: number of substrate meshpoints in x is xSteps - 1
        if y % 2 == 0:
            for x in range(xSteps - 1):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0
        # if y is odd: number of substrate meshpoints in x is xSteps
        else:
            for x in range(xSteps):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

    return
