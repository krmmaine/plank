from numpy import linspace
from matplotlib import pyplot
from mpl_toolkits import mplot3d            # Still need this even though is seems unnecessary


def createGraph(ySubstrate, xSteps, vegf, xVector, yVector):

    zVector = []

    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                zVector.append(vegf[y][x])
        else:
            for x in range(xSteps):
                zVector.append(vegf[y][x])

    ax = pyplot.axes(projection='3d')
    ax.scatter(xVector, yVector, zVector, c=zVector, cmap='viridis', linewidth=0.5)
    ax = pyplot.axes(projection='3d')
    ax.plot_trisurf(xVector, yVector, zVector,
                    cmap='viridis', edgecolor='none')
    pyplot.show()

    return
