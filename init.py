def init(xSteps, totNumCells, xPos, yPos, occupied, deathTime, numTimeSteps):
    cellNum = 0
    # seed initial cells
    for x in range(int(0.5*xSteps/totNumCells), xSteps, int(xSteps/totNumCells)):
        xPos[cellNum][0] = x
        yPos[cellNum][0] = 0
        occupied[0][x] += 1
        deathTime[cellNum] = numTimeSteps - 1
        cellNum += 1

    return
