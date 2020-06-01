from numpy import zeros


def updateFib(ySubstrate, xSteps, densityScale, occupiedOld, fib, fibOld, k, pro, tolerance, h):
    f = zeros((ySubstrate, xSteps))
    k4 = 154.3209877
    k5 = 1740294.403
    k6 = 0.012048193
    relax2 = 1
    k22 = 0.0001
    k23 = 154.3209877

    # capillary
    for x in range(xSteps-1):
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x+1]) / 2
        fibOld[0][x] = fib[0][x]
        fib[0][x] = fib[0][x] + k * (k4 * fib[0][x] * (1 - fib[0][x]) * density - k5 * pro[0][x] * fib[0][x] /
                                     (1 + k6 * fib[0][x]))    # EQ 48 amd EQ 51
        if fib[0][x] < 0:
            fib[0][x] = 0
        if fib[0][x] > 1:
            fib[0][x] = 1

    # initialize f: guess
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                f[y][x] = fib[y][x]
        else:
            for x in range(xSteps):
                f[y][x] = fib[y][x]

    intol = 0
    while intol == 0:
        intol = 1

    # boundary at y max
        fOld = f[ySubstrate - 1][0]
        f[ySubstrate - 1][0] = f[ySubstrate - 1][1]
        if f[ySubstrate - 1][0] - fOld > tolerance or f[ySubstrate - 1][0] - fOld < -tolerance:
            intol = 0

        for x in range(1, xSteps-2, 1):
            fOld = f[ySubstrate - 1][x]
            f[ySubstrate - 1][x] = f[ySubstrate - 3][x]
            if f[ySubstrate - 1][x] - fOld > tolerance or f[ySubstrate - 1][x] - fOld < -tolerance:
                intol=0

        fOld = f[ySubstrate - 1][xSteps - 2]
        f[ySubstrate - 1][xSteps - 2] = f[ySubstrate - 1][xSteps - 3]

        if f[ySubstrate - 1][xSteps - 2] - fOld > tolerance or f[ySubstrate - 1][xSteps - 2] - fOld < -tolerance:
            intol = 0

    # boundary at y max - 1
        fOld = f[ySubstrate-2][0]
        f[ySubstrate-2][0] = f[ySubstrate-2][1]
        if f[ySubstrate-2][0] - fOld > tolerance or f[ySubstrate-2][0] - fOld < -tolerance:
            intol=0

        for x in range(1, xSteps-1, 1):
            fOld = f[ySubstrate-2][x]
            f[ySubstrate-2][x] = f[ySubstrate-4][x]
            if f[ySubstrate-2][x] - fOld > tolerance or f[ySubstrate-2][x] - fOld < -tolerance:
                intol=0

        fOld = f[ySubstrate-2][xSteps-1]
        f[ySubstrate-2][xSteps-1] = f[ySubstrate-2][xSteps-2]

        if f[ySubstrate-2][xSteps-1] - fOld > tolerance or f[ySubstrate-2][xSteps-1] - fOld < -tolerance:
            intol = 0

    # cycle through interior rows
        for y in range(ySubstrate-3, 2, -1):
            fOld = f[y][0]
            f[y][0] = f[y][1]

            if f[y][0] - fOld > tolerance or f[y][0] - fOld < -tolerance:
                intol = 0

            if y % 2 == 0:
                for x in range(1, xSteps-2, 1):
                    fOld = f[y][x]
                    f[y][x] = relax2 * k /(h*h+2*k22*k) * (0.5*k22*
                                                           (f[y][x+1]+f[y][x-1]+f[y+2][x]+f[y-2][x]
                                                            +fib[y][x+1]+fib[y][x-1]+fib[y+2][x]+fib[y-2][x])
                                                           + (h*h/k - 2*k22 +
                                                              k23*h*h*(1-fib[y][x]) - k5*h*h*pro[y][x]/(1+k6*fib[y][x]))
                                                           * fib[y][x]) + (1-relax2)*f[y][x]    # EQ 55 and EQ 59 pg 179
                    if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                        intol = 0

                fOld = f[y][xSteps-2]
                f[y][xSteps-2] = f[y][xSteps-3]
                if f[y][xSteps-2] - fOld > tolerance or f[y][xSteps-2] - fOld < -tolerance:
                    intol = 0

            else:
                for x in range(1, xSteps-1, 1):
                    fOld = f[y][x]
                    f[y][x] = relax2 * k / (h * h + 2 * k22 * k) * (0.5 * k22 *
                                                                (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                                                + fib[y][x+1] + fib[y][x-1] + fib[y+2][x] + fib[y-2][x])
                                                                + (h * h / k - 2 * k22 +
                                                                    k23 * h * h * (1 - fib[y][x]) - k5 * h * h *
                                                                    pro[y][x] / (1 + k6 * fib[y][x]))* fib[y][x]) \
                                                                    + (1 - relax2) * f[y][x]  # EQ 55 and EQ 59 pg 179

                    if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                        intol=0
                fOld = f[y][xSteps-1]
                f[y][xSteps-1] = f[y][xSteps-2]
                if f[y][xSteps-1] - fOld > tolerance or f[y][xSteps-1] - fOld < -tolerance:
                    intol = 0
                    
    # boundary row at y = 2
        fOld = f[2][0]
        f[2][0] = f[2][1]
        if f[2][0] - fOld > tolerance or f[2][0] - fOld < -tolerance:
            intol=0
        for x in range(1, xSteps-2, 1):
            fOld = f[2][x]
            f[2][x] = f[4][x]
            if f[2][x] - fOld > tolerance or f[2][x] - fOld < -tolerance:
                intol=0
        fOld = f[2][xSteps-2]
        f[2][xSteps-2] = f[2][xSteps-3]
        if f[2][xSteps-2] - fOld > tolerance or f[2][xSteps-2] - fOld < -tolerance:
            intol = 0

    # boundary row at y = 1
        fOld = f[1][0]
        f[1][0] = f[1][1]
        if f[1][0] - fOld > tolerance or f[1][0] - fOld < -tolerance:
            intol=0

        for x in range(1, xSteps-1, 1):
            fOld = f[1][x]
            f[1][x] = f[3][x]
            if f[1][x] - fOld > tolerance or f[1][x] - fOld < -tolerance:
                intol=0
        fOld = f[1][xSteps-1]
        f[1][xSteps-1] = f[1][xSteps-2]

        if f[1][xSteps-1] - fOld > tolerance or f[1][xSteps-1] - fOld < -tolerance:
            intol = 0

# Cycle through substrate meshpoints and set fibronectin at time step j+1
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                fibOld[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1
        else:
            for x in range(xSteps):
                fibOld[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1

    return
