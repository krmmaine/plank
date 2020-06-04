

def updatePro(ySubstrate, xSteps, densityScale, occupiedOld, pro, proOld, k, vegfOld):
    k1 = 51923.07639
    k3 = 3166.6666667
    totProCap = 0

# update protease concentration in capillary
    for x in range(xSteps-1):
        # Average density at cell meshpoints to the right and left of substrate meshpoint at time j
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x+1]) / 2
        proOld[0][x] = pro[0][x]
        # Use equation 47 to update protease concentration in the capillary
        # use vegfOld because vegf has already been updated this time step
        pro[0][x] = pro[0][x] + k * (k1 * vegfOld[0][x] * density / (vegfOld[0][x]+1) - k3 * pro[0][x])
        if pro[0][x] < 0:
            pro[0][x] = 0
        totProCap = totProCap + pro[0][x]
    totProCap = totProCap / xSteps

# update protease concentratino at boundary y = 1
    for x in range(xSteps):
        # don't take into account cells still in the capillary
        density = densityScale ** 2 * occupiedOld[1][x]
        proOld[1][x] = pro[1][x]
        # use equation 54 (same as 47 but for the ECM instead of the capillary)
        pro[1][x] = pro[1][x] + k * (k1 * vegfOld[1][x] * density / (vegfOld[1][x]+1) - k3 * pro[1][x])
        if pro[1][x] < 0:
            pro[1][x] = 0
            
# cycle through substrate meshpoints
    for y in range(2, ySubstrate, 1):
        # If y is even, number of substrate meshpoints in x-direction is xSteps-1
        if y % 2 == 0:        
            for x in range(xSteps-1):
                # Average density at cell meshpoints to the right and left of substrate meshpoint at time j
                density = densityScale ** 2 * (occupiedOld[y//2][x] + occupiedOld[y//2][x+1]) / 2
                proOld[y][x] = pro[y][x]
                # use equation 54 (same as 47 but for the ECM instead of the capillary)
                pro[y][x] = pro[y][x] + k * (k1 * vegfOld[y][x] * density / (vegfOld[y][x]+1) - k3 * pro[y][x])
                if pro[y][x] < 0:
                    pro[y][x] = 0
        # If y is odd, number of substrate meshpoints in x-direction is xSteps
        else:
            for x in range(xSteps):
                # Average density at cell meshpoints above and below substrate meshpoint at time j
                density = densityScale ** 2 * (occupiedOld[(y-1)//2][x] + occupiedOld[(y+1)//2][x]) / 2
                proOld[y][x] = pro[y][x]
                # use equation 54 (same as 47 but for the ECM instead of the capillary)
                pro[y][x] = pro[y][x] + k * (k1 * vegfOld[y][x] * density / (vegfOld[y][x]+1) - k3 * pro[y][x])
                if pro[y][x] < 0:
                    pro[y][x] = 0

    return
