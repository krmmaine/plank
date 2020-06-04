def lamda(xLength, stepSize):
    # non-dimensionalized equation on page 150
    # can't use lambda because that is a python keyword
    lamda = (xLength**2) / (stepSize**2)
    return lamda
