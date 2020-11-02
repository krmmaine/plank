from numpy import linspace
from matplotlib import pyplot
from mpl_toolkits import mplot3d            # Still need this even though is seems unnecessary
from matplotlib import cm


def createGraph(ySubstrate, xSteps, vegf, fib, pro, xVector, yVector, workspace, time):

    VEGFzVector = []

    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                VEGFzVector.append(vegf[y][x])
        else:
            for x in range(xSteps):
                VEGFzVector.append(vegf[y][x])

    FIBROzVector = []

    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                FIBROzVector.append(fib[y][x])
        else:
            for x in range(xSteps):
                FIBROzVector.append(fib[y][x])

    PROzVector = []

    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                PROzVector.append(pro[y][x])
        else:
            for x in range(xSteps):
                PROzVector.append(pro[y][x])

    fileName = "Time = " + str(time) + ".pdf"
    file = open(fileName, "a+")

    fig = pyplot.figure(figsize=pyplot.figaspect(2.7))

    ax = fig.add_subplot(4, 1, 1)
    ax.imshow(workspace)
    cm.get_cmap("jet")

    ax = fig.add_subplot(4, 1, 2, projection='3d')
    ax.plot_trisurf(xVector, yVector, VEGFzVector, cmap='viridis', edgecolor='none')

    ax = fig.add_subplot(4, 1, 3, projection='3d')
    ax.plot_trisurf(xVector, yVector, FIBROzVector, cmap='viridis', edgecolor='none')

    ax = fig.add_subplot(4, 1, 4, projection='3d')
    ax.plot_trisurf(xVector, yVector, PROzVector, cmap='viridis', edgecolor='none')

    # 2D VEGF graph
    # ax = fig.add_subplot(5, 1, 5)
    # ax.imshow(vegf)
    # cm.get_cmap("jet")

    pyplot.savefig(fileName)
    file.close()

    return
