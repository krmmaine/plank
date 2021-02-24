import numpy as np
from math import cos
from math import pi

from xCoordinate import xCoordinate


def newUpdateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength, v0, Dv):

    sigma = 1.514705513 * (10 ** -3)
    v1 = 0.007692308

    k35 = (v0 * sigma * v1) / Dv

    Dp = 3.6 * (10 ** -6)

    k21 = Dv / Dp

    k1 = 51923.07639
    k2 = 694.4444444
    m0 = 12
    relax = 1.45
    k33 = 0.00036

    # update VEGF concentration in capillary
    vegfOld[0, :] = vegf[0, :]
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
        # update VEGF concentration in capillary using equation 46
        vegf[0][x] = vegf[0][x] + k * (-k1 * vegf[0][x] * density / (1 + vegf[0][x]) + k2 * VEGFdiff)
        if vegf[0][x] < 0:
            vegf[0][x] = 0

    v = np.copy(vegf)

    # Keep iterating until the value of v at each meshpoint changes by less than a certain tolerance
    intol = 0
    while intol == 0:
        intol = 1

        for y in range(ySubstrate - 1, 0, -1):

            if y % 2 == 0:
                stop = xSteps - 2
            else:
                stop = xSteps - 1

            # update VEGF concentration at boundary at maximum y. and y-1 includes source
            if y == ySubstrate - 1 or y == ySubstrate - 2:
                # using equation 70 derivation on page 179 calculate vegf at x = 0
                vOld = v[y][0]
                v[y][0] = v[y][1]
                # if the change in any v and vOld is greater than the tolerance then this loop will continue to run
                if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                    intol = 0

                for x in range(1, stop, 1):
                    vOld = v[y][x]
                    # use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary
                    v[y][x] = k35 * h * ((1 - cos(2 * pi * xCoordinate(x, y, xSteps, xLength))) ** m0) + v[y - 2][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol = 0

                # using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][stop]
                v[y][stop] = v[y][stop - 1]
                if v[y][stop] - vOld > tolerance or v[y][stop] - vOld < -tolerance:
                    intol = 0

            # update VEGF concentration at capillary boundary
            elif y == 1 or y == 2:
                # using equation 70 derivation on page 179 calculate vegf at x = 0
                vOld = v[y][0]
                v[y][0] = v[y][1]
                if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                    intol = 0

                for x in range(1, stop, 1):
                    if y == 2:
                        average = vegfOld[0][x]
                    else:
                        average = (vegfOld[0][x - 1] + vegfOld[0][x] / 2)
                    vOld = v[y][x]
                    # use equation 61 see derivation on page 179 of paper and pg 137 of notes
                    # use vegfOld because capillary values have already been updated
                    v[y][x] = 1 / (k33 + h) * (k33 * v[y + 2][x] + h * average)
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        intol = 0

                # using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][stop]
                v[y][stop] = v[y][stop - 1]
                if v[y][stop] - vOld > tolerance or v[y][stop] - vOld < -tolerance:
                    intol = 0

            else:
                # using equation 70 derivation on page 179 calculate vegf at x = 0
                vOld = v[y][0]
                v[y][0] = v[y][1]
                if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                    intol = 0

                for x in range(1, stop, 1):
                    # average density of cell meshpoint to the right and left
                    if y % 2 == 0:
                        density = densityScale * (ySubstrate / 2 - 1) * (
                                occupiedOld[y // 2][x] + occupiedOld[y // 2][x + 1]) / 2
                    else:
                        # average density of cell meshpoints above and below
                        # y // 2 because there are twice as many points in the y direction because substrate meshpoints are at 1/2
                        density = densityScale * (ySubstrate / 2 - 1) * (
                                occupiedOld[(y - 1) // 2][x] + occupiedOld[(y + 1) // 2][x]) / 2
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
                vOld = v[y][stop]
                v[y][stop] = v[y][stop-1]
                if v[y][stop] - vOld > tolerance or v[y][stop] - vOld < -tolerance:
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

    return vegfOld, vegf
