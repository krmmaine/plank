from tau import tau


def pMove(x, y, direction, pro, fib, vegf, xSteps, ySteps, lamda, k, movement):
    T = [0, 0, 0, 0]

    # LEFT
    if x > 0:
        if y == 0:
            T[0] = tau(pro[2*y][x-1], fib[2*y][x-1], vegf[2*y][x-1], y)
        elif movement[len(movement)-1] == 2:
            T[0] = 0
        else:
            T[0] = tau(pro[2 * y][x - 1], fib[2 * y][x - 1], vegf[2 * y][x - 1], y)
    else:
        T[0] = 0

    # RIGHT
    if x < xSteps-1:
        if y == 0:
            T[1] = tau(pro[2*y][x], fib[2*y][x], vegf[2*y][x], y)       # don't need to add 1 to x because the x substrate should be offset by a half step from the cell matrix meshpoints
        elif movement[len(movement)-1] == 3:
            T[1] = 0
        else:
            T[1] = tau(pro[2 * y][x], fib[2 * y][x], vegf[2 * y][x], y)
    else:
        T[1] = 0

    # UP
    if y > 1 and movement[len(movement)-1] != 0:
        T[2] = tau(pro[2*y-1][x], fib[2*y-1][x], vegf[2*y-1][x], y)
    else:
        T[2] = 0

    # DOWN
    if 0 < y < ySteps - 1 and movement[len(movement)-1] != 1:
        T[3] = tau(pro[2*y+1][x], fib[2*y+1][x], vegf[2*y+1][x], y)
    else:
        T[3] = 0

    # probability equations from page 153 and 152
    if y == 0:
        if T[0]+T[1] != 0:
            pmove = 2 * lamda * k * T[direction] / (T[0]+T[1])
        else:
            pmove = 0
    else:
        if T[0]+T[1]+T[2]+T[3] != 0:
            pmove = 4 * lamda * k * T[direction] / (T[0]+T[1]+T[2]+T[3])
        else:
            pmove = 0

    return pmove, T
