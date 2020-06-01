def xCoordinate(x, y, xSteps, xLength):
    if y % 2 == 0:
        return ((x + 0.5) / xSteps) * xLength
    else:
        return (x / xSteps) * xLength