from numpy import zeros
from math import cos
from math import pi

from xCoordinate import xCoordinate

def updateSubstrates(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, fib, fibOld, pro, proOld, k, tolerance, h, xLength):
    k1 = 51923.07639
    k2 = 694.4444444
    k3 = 3166.6666667
    k4 = 154.3209877
    k5 = 1740294.403
    k6 = 0.012048193
    k21 = 10
    k22 = 0.0001
    k23 = 154.3209877
    k33 = 0.00036
    k35 = 0.0129462015
    m0 = 12
    relax = 1.45
    relax2 = 1
    totProCap = 0

    v = zeros((ySubstrate, xSteps))
    f = zeros((ySubstrate, xSteps))

    #Update VEGF, fib, and pro concentrations in the capillary
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
        fibOld[0][x] = fib[0][x]
        # update VEGF concentration in capillary using equation 46
        vegf[0][x] = vegf[0][x] + k * (-k1 * vegf[0][x] * density / (1 + vegf[0][x]) + k2 * VEGFdiff)
        # calculate fibronectin concentration in capillary using equation 48 and equation 51
        fib[0][x] = fib[0][x] + k * (k4 * fib[0][x] * (1 - fib[0][x]) * density - k5 * pro[0][x] * fib[0][x] /
                                     (1 + k6 * fib[0][x]))
        if vegf[0][x] < 0:
            vegf[0][x] = 0
        if fib[0][x] < 0:
            fib[0][x] = 0
        # one is the initial and max value of the fibronectin concentration
        if fib[0][x] > 1:
            fib[0][x] = 1

        # Average density at cell meshpoints to the right and left of substrate meshpoint at time j
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x + 1]) / 2
        proOld[0][x] = pro[0][x]
        # Use equation 47 to update protease concentration in the capillary
        # use vegfOld because vegf has already been updated this time step
        pro[0][x] = pro[0][x] + k * (k1 * vegfOld[0][x] * density / (vegfOld[0][x] + 1) - k3 * pro[0][x])
        if pro[0][x] < 0:
            pro[0][x] = 0
        totProCap = totProCap + pro[0][x]
    totProCap = totProCap / xSteps

    # initialize v and f: value of the previous time step of vegf and fib, respectively
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                v[y][x] = vegf[y][x]
            for x in range(xSteps - 1):
                f[y][x] = fib[y][x]
        else:
            for x in range(xSteps):
                v[y][x] = vegf[y][x]
            for x in range(xSteps):
                f[y][x] = fib[y][x]

# Keep iterating until the value of v at each meshpoint changes by less than a certain tolerance
    vintol = 0 # tolerance for vegf
    fintol = 0 # tolerance for fib
    # when intol is 0, a check has failed and the variable is not yet within tolerance.
    # when intol is 1, the variable is temporarily within tolerance, but may not pass a check later in the loop
    # when intol is 2, the variable has passed all checks and no longer needs to be updated in the loop
    count = 0
    while vintol != 2 or fintol != 2:
        if vintol != 2:
            vintol = 1
        if fintol != 2:
            fintol = 1
        #print(count, 0, vintol)

        # update VEGF concentration at boundary at maximum y. includes source
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        if vintol != 2:
            vOld = v[ySubstrate - 1][0]
            v[ySubstrate - 1][0] = v[ySubstrate - 1][1]
            # if the change in any v and vOld is greater than the tolerance then this loop will continue to run
            if v[ySubstrate - 1][0] - vOld > tolerance or v[ySubstrate - 1][0] - vOld < -tolerance:
                vintol = 0
        #print(count, 1, vintol)

        # update fibronectin concentration at boundary at maximum y. includes source
        # using equation 71 derivation on page 179 calculate fib at x = 0
        if fintol != 2:
            fOld = f[ySubstrate - 1][0]
            f[ySubstrate - 1][0] = f[ySubstrate - 1][1]
            # if the change in any f and fOld is ever greater than the tolerance then this loop will continue to run
            if f[ySubstrate - 1][0] - fOld > tolerance or f[ySubstrate - 1][0] - fOld < -tolerance:
                fintol = 0

        for x in range(1, xSteps - 2, 1):
            if vintol != 2:
                vOld = v[ySubstrate - 1][x]
                # use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary
                v[ySubstrate - 1][x] = k35 * h * ((1 - cos(2 * pi * xCoordinate(x, ySubstrate - 1, xSteps, xLength))) ** m0) \
                                       + v[ySubstrate - 3][x]
                if v[ySubstrate - 1][x] - vOld > tolerance or v[ySubstrate - 1][x] - vOld < -tolerance:
                    vintol = 0

            if fintol != 2:
                fOld = f[ySubstrate - 1][x]
                # use EQ 66 and derivation on page 179 to update fib concentration at upper boundary
                f[ySubstrate - 1][x] = f[ySubstrate - 3][x]
                if f[ySubstrate - 1][x] - fOld > tolerance or f[ySubstrate - 1][x] - fOld < -tolerance:
                    fintol = 0
        #print(count, 2, vintol)

        # using equation 70 derivation on page 179 calculate vegf at x = max
        if vintol != 2:
            vOld = v[ySubstrate - 1][xSteps - 2]
            v[ySubstrate - 1][xSteps - 2] = v[ySubstrate - 1][
                xSteps - 3]  # minus 3 because row is even
            if v[ySubstrate - 1][xSteps - 2] - vOld > tolerance or v[ySubstrate - 1][xSteps - 2] - vOld < -tolerance:
                vintol = 0
        #print(count, 3, vintol)

        # using equation 71 derivation on page 179 calculate fib at x = max
        if fintol != 2:
            fOld = f[ySubstrate - 1][xSteps - 2]
            f[ySubstrate - 1][xSteps - 2] = f[ySubstrate - 1][xSteps - 3]
            if f[ySubstrate - 1][xSteps - 2] - fOld > tolerance or f[ySubstrate - 1][xSteps - 2] - fOld < -tolerance:
                fintol = 0

        #update VEGF concentration at boundary at maximum y - 1
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        if vintol != 2:
            vOld = v[ySubstrate - 2][0]
            v[ySubstrate - 2][0] = v[ySubstrate - 2][1]
            if v[ySubstrate - 2][0] - vOld > tolerance or v[ySubstrate - 2][0] - vOld < -tolerance:
                vintol = 0
        #print(count, 4, vintol)

        # update fib concentration at boundary at maximum y - 1
        # using equation 71 derivation on page 179 calculate fib at x = 0
        if fintol != 2:
            fOld = f[ySubstrate - 2][0]
            f[ySubstrate - 2][0] = f[ySubstrate - 2][1]
            if f[ySubstrate - 2][0] - fOld > tolerance or f[ySubstrate - 2][0] - fOld < -tolerance:
                fintol = 0

        for x in range(1, xSteps - 1, 1):
            if vintol != 2:
                vOld = v[ySubstrate - 2][x]
                # use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary - 1
                v[ySubstrate - 2][x] = k35 * h * ((1 - cos(2 * pi * xCoordinate(x, ySubstrate - 2, xSteps, xLength))) ** m0) \
                                       + v[ySubstrate - 4][x]  # -4 because substrate points are at half meshpoints
                if v[ySubstrate - 2][x] - vOld > tolerance or v[ySubstrate - 2][x] - vOld < -tolerance:
                    vintol = 0

            if fintol != 2:
                fOld = f[ySubstrate - 2][x]
                # use EQ 66 and derivation on page 179 to update fib concentration at upper boundary - 1
                f[ySubstrate - 2][x] = f[ySubstrate - 4][x]
                if f[ySubstrate - 2][x] - fOld > tolerance or f[ySubstrate - 2][x] - fOld < -tolerance:
                    fintol = 0
        #print(count, 5, vintol)

        # using equation 70 derivation on page 179 calculate vegf at x = max
        if vintol != 2:
            vOld = v[ySubstrate - 2][xSteps - 1]
            v[ySubstrate - 2][xSteps - 1] = v[ySubstrate - 2][xSteps - 2]
            if v[ySubstrate - 2][xSteps - 1] - vOld > tolerance or v[ySubstrate - 2][xSteps - 1] - vOld < -tolerance:
                vintol = 0
        #print(count, 6, vintol)

        # using equation 71 derivation on page 179 calculate fib at x = max
        if fintol != 2:
            fOld = f[ySubstrate - 2][xSteps - 1]
            f[ySubstrate - 2][xSteps - 1] = f[ySubstrate - 2][xSteps - 2]
            if f[ySubstrate - 2][xSteps - 1] - fOld > tolerance or f[ySubstrate - 2][xSteps - 1] - fOld < -tolerance:
                fintol = 0

    # Cycle through interior rows and calculate new VEGF concentration
        for y in range(ySubstrate - 3, 2, -1):
            # using equation 70 derivation on page 179 calculate vegf at x = 0
            if vintol != 2:
                vOld = v[y][0]
                v[y][0] = v[y][1]
                if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                    vintol = 0
            #print(count, 7, vintol)

            # using equation 71 derivation on page 179 calculate fib at x = 0
            if fintol != 2:
                fOld = f[y][0]
                f[y][0] = f[y][1]
                if f[y][0] - fOld > tolerance or f[y][0] - fOld < -tolerance:
                    fintol = 0

            # if row is even number of substrate meshpoints in x is xsteps-1
            if y % 2 == 0:
                for x in range(1, xSteps - 2, 1):  # minus 2 because last meshpoint has a boundary condition
                    # densityScale is squared because density in ECM is per unit area not length
                    # average density of cell meshpoint to the right and left
                    density = densityScale * (ySubstrate/2-1) * (occupiedOld[y // 2][x] + occupiedOld[y // 2][x + 1]) / 2
                    if vintol != 2:
                        vOld = v[y][x]
                        # Approximate equation 53 using the crank-nicolson method see derivation on page 178
                        # relax is the successive over-relaxation term
                        v[y][x] = relax / (h * h + 2 * k21 * k) * (
                                    0.5 * k21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x] +
                                                     vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                    + (h * h - 2 * k21 * k - h * h * k * k1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                                  + (1 - relax) * v[y][x]
                        if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                            vintol = 0

                    if fintol != 2:
                        fOld = f[y][x]
                        # Approximate equation 55 and equation 59 using the crank-nicolson method see derivation on page 179
                        # relax2 is the successive over-relaxation term
                        f[y][x] = relax2 * k / (h * h + 2 * k22 * k) * (0.5 * k22 *
                                                                        (f[y][x + 1] + f[y][x - 1] + f[y + 2][x] + f[y - 2][
                                                                            x]
                                                                         + fib[y][x + 1] + fib[y][x - 1] + fib[y + 2][x] +
                                                                         fib[y - 2][x])
                                                                        + (h * h / k - 2 * k22 +
                                                                           k23 * h * h * (1 - fib[y][x]) - k5 * h * h *
                                                                           pro[y][x] / (1 + k6 * fib[y][x]))
                                                                        * fib[y][x]) + (1 - relax2) * f[y][x]
                        if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                            fintol = 0
                #print(count, 8, vintol)

                # using equation 70 derivation on page 179 calculate vegf at x = max
                if vintol != 2:
                    vOld = v[y][xSteps - 2]
                    v[y][xSteps - 2] = v[y][xSteps - 3]
                    if v[y][xSteps - 2] - vOld > tolerance or v[y][xSteps - 2] - vOld < -tolerance:
                        vintol = 0
                #print(count, 9, vintol)

                # using equation 71 derivation on page 179 calculate fib at x = max
                if fintol != 2:
                    fOld = f[y][xSteps - 2]
                    f[y][xSteps - 2] = f[y][xSteps - 3]
                    if f[y][xSteps - 2] - fOld > tolerance or f[y][xSteps - 2] - fOld < -tolerance:
                        fintol = 0

            else:
                # if row is odd number of substrate meshpoints in x is nn
                for x in range(1, xSteps - 1, 1):
                    # average density of cell meshpoints above and below
                    # y // 2 because there are twice as many points in the y direction because substrate meshpoints are at 1/2
                    density = densityScale * (ySubstrate/2-1) * (occupiedOld[(y - 1) // 2][x] + occupiedOld[(y + 1) // 2][x]) / 2
                    if vintol != 2:
                        vOld = v[y][x]
                        v[y][x] = relax / (h * h + 2 * k21 * k) * (
                                    0.5 * k21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x]
                                                     + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                    + (h * h - 2 * k21 * k - h * h * k * k1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                                  + (1 - relax) * v[y][x]
                        if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                            vintol = 0

                    if fintol != 2:
                        fOld = f[y][x]
                        # Approximate equation 55 and equation 59 using the crank-nicolson method see derivation on page 179
                        # relax2 is the successive over-relaxation term
                        f[y][x] = relax2 * k / (h * h + 2 * k22 * k) * (0.5 * k22 *
                                                                        (f[y][x + 1] + f[y][x - 1] + f[y + 2][x] + f[y - 2][
                                                                            x]
                                                                         + fib[y][x + 1] + fib[y][x - 1] + fib[y + 2][x] +
                                                                         fib[y - 2][x])
                                                                        + (h * h / k - 2 * k22 +
                                                                           k23 * h * h * (1 - fib[y][x]) - k5 * h * h *
                                                                           pro[y][x] / (1 + k6 * fib[y][x])) * fib[y][x]) \
                                  + (1 - relax2) * f[y][x]
                        if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                            fintol = 0
                #print(count, 10, vintol)

                # using equation 70 derivation on page 179 calculate vegf at x = max
                if vintol != 2:
                    vOld = v[y][xSteps - 1]
                    v[y][xSteps - 1] = v[y][xSteps - 2]
                    if v[y][xSteps - 1] - vOld > tolerance or v[y][xSteps - 1] - vOld < -tolerance:
                        vintol = 0
                #print(count, 11, vintol)

                # using equation 71 derivation on page 179 calculate fib at x = max
                if fintol != 2:
                    fOld = f[y][xSteps - 1]
                    f[y][xSteps - 1] = f[y][xSteps - 2]
                    if f[y][xSteps - 1] - fOld > tolerance or f[y][xSteps - 1] - fOld < -tolerance:
                        fintol = 0

        # calculate VEGF concentration at boundary row y = 2
        # using equation 70 derivation on page 179 calculate vegf at x = 0
        if vintol != 2:
            vOld = v[2][0]
            v[2][0] = v[2][1]
            if v[2][0] - vOld > tolerance or v[2][0] - vOld < -tolerance:
                vintol = 0
        #print(count, 12, vintol)

        # calculate fibronectin concentration at boundary row y = 2
        # using equation 71 derivation on page 179 calculate fib at x = 0
        if fintol != 2:
            fOld = f[2][0]
            f[2][0] = f[2][1]
            if f[2][0] - fOld > tolerance or f[2][0] - fOld < -tolerance:
                fintol = 0

        for x in range(1, xSteps - 2, 1):
            if vintol != 2:
                vOld = v[2][x]
                # use equation 61 see derivation on page 179 of paper and pg 137 of notes
                # use vegfOld because capillary values have already been updated
                v[2][x] = 1 / (k33 + h) * (k33 * v[4][x] + h * vegfOld[0][x])
                if v[2][x] - vOld > tolerance or v[2][x] - vOld < -tolerance:
                    vintol = 0

            if fintol != 2:
                fOld = f[2][x]
                # use EQ 62 and derivation on page 179 to update fib concentration at lower boundary
                f[2][x] = f[4][x]
                if f[2][x] - fOld > tolerance or f[2][x] - fOld < -tolerance:
                    fintol = 0
        #print(count, 13, vintol)

        # using equation 70 derivation on page 179 calculate vegf at x = max
        if vintol != 2:
            vOld = v[2][xSteps - 2]
            v[2][xSteps - 2] = v[2][xSteps - 3]
            if v[2][xSteps - 2] - vOld > tolerance or v[2][xSteps - 2] - vOld < -tolerance:
                vintol = 0
        #print(count, 14, vintol)

        # using equation 71 derivation on page 179 calculate fib at x = max
        if fintol != 2:
            fOld = f[2][xSteps - 2]
            f[2][xSteps - 2] = f[2][xSteps - 3]
            if f[2][xSteps - 2] - fOld > tolerance or f[2][xSteps - 2] - fOld < -tolerance:
                fintol = 0

        # calculate VEGF concentratin at boundary row y = 1: capillary wall
        if vintol != 2:
            vOld = v[1][0]
            v[1][0] = v[1][1]
            if v[1][0] - vOld > tolerance or v[1][0] - vOld < -tolerance:
                vintol = 0
        #print(count, 15, vintol)

        # using equation 71 derivation on page 179 calculate fib at x = 0
        if fintol != 2:
            fOld = f[1][0]
            f[1][0] = f[1][1]
            if f[1][0] - fOld > tolerance or f[1][0] - fOld < -tolerance:
                fintol = 0

        for x in range(1, xSteps - 1, 1):
            if vintol != 2:
                vOld = v[1][x]
                # Take average because there is no meshpoint in the capillary directly below
                v[1][x] = 1 / (k33 + h) * (k33 * v[3][x] + h * (vegfOld[0][x - 1] + vegfOld[0][x]) / 2)
                if v[1][x] - vOld > tolerance or v[1][x] - vOld < -tolerance:
                    vintol = 0

            if fintol != 2:
                # use EQ 62 and derivation on page 179 to update fib concentration at lower boundary
                fOld = f[1][x]
                f[1][x] = f[3][x]
                if f[1][x] - fOld > tolerance or f[1][x] - fOld < -tolerance:
                    fintol = 0
        #print(count, 16, vintol)

        # using equation 70 derivation on page 179 calculate vegf at x = max
        if vintol != 2:
            vOld = v[1][xSteps - 1]
            v[1][xSteps - 1] = v[1][xSteps - 2]
            if v[1][xSteps - 1] - vOld > tolerance or v[1][xSteps - 1] - vOld < -tolerance:
                vintol = 0
        #print(count, 17, vintol)

        if fintol != 2:
            # using equation 71 derivation on page 179 calculate fib at x = max
            fOld = f[1][xSteps - 1]
            f[1][xSteps - 1] = f[1][xSteps - 2]
            if f[1][xSteps - 1] - fOld > tolerance or f[1][xSteps - 1] - fOld < -tolerance:
                fintol = 0

        if vintol == 1:
            vintol = 2
        if fintol == 1:
            fintol = 2
        #print(count, 18, vintol)
        count += 1

    # Cycle through substrate meshpoints and set VEGF and fib at time step j+1
    for y in range(1, ySubstrate, 1):
        # if y is even: number of substrate meshpoints in x is xSteps - 1
        if y % 2 == 0:
            for x in range(xSteps - 1):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

                fibOld[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1

        # if y is odd: number of substrate meshpoints in x is xSteps
        else:
            for x in range(xSteps):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

                fibOld[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1

    return