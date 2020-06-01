

def updatePro(ySubstrate, xSteps, densityScale, occupiedOld, pro, proOld, k, vegfOld):
    k1 = 51923.07639
    k3 = 3166.6666667
    totProCap = 0

# capillary at y = 0
    for x in range(xSteps-1):
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x+1]) / 2
        proOld[0][x] = pro[0][x]
        pro[0][x] = pro[0][x] + k * (k1 * vegfOld[0][x] * density / (vegfOld[0][x]+1) - k3 * pro[0][x])     # EQ 47
        if pro[0][x] < 0:
            pro[0][x] = 0
        totProCap = totProCap + pro[0][x]
    totProCap = totProCap / xSteps

# boundary at y = 1
    for x in range(xSteps):
        density = densityScale ** 2 * occupiedOld[1][x]
        proOld[1][x] = pro[1][x]
        pro[1][x] = pro[1][x] + k * (k1 * vegfOld[1][x] * density / (vegfOld[1][x]+1) - k3 * pro[1][x])     # EQ 47
        if pro[1][x] < 0:
            pro[1][x] = 0
            
# cycle through substrate meshpoints
    for y in range(2, ySubstrate, 1):        
        if y % 2 == 0:        
            for x in range(xSteps-1):        
                density = densityScale ** 2 * (occupiedOld[y//2][x] + occupiedOld[y//2][x+1]) / 2
                proOld[y][x] = pro[y][x]
                pro[y][x] = pro[y][x] + k * (k1 * vegfOld[y][x] * density / (vegfOld[y][x]+1) - k3 * pro[y][x]) # EQ 47
                if pro[y][x] < 0:
                    pro[y][x] = 0
        
        else:
            for x in range(xSteps):
                density = densityScale ** 2 * (occupiedOld[(y-1)//2][x] + occupiedOld[(y+1)//2][x]) / 2
                proOld[y][x] = pro[y][x]
                pro[y][x] = pro[y][x] + k * (k1 * vegfOld[y][x] * density / (vegfOld[y][x]+1) - k3 * pro[y][x]) # EQ 47
                if pro[y][x] < 0:
                    pro[y][x] = 0

    return
