def move(cell, time, stay, left, right, down, random, yPos, xPos, occupied, fib, vegf, pro, movedUp, movedDown, movedLeft, movedRight):

    # UP
    if random > stay + left + right + down:
        yPos[cell][time+1] = yPos[cell][time] + 1
        print("moved up")
        movedUp += 1
    # DOWN
    elif random > stay + left + right:
        yPos[cell][time+1] = yPos[cell][time] - 1
        print("moved down")
        movedDown += 1
    # RIGHT
    elif random > stay + left:
        xPos[cell][time+1] = xPos[cell][time] + 1
        print("moved right")
        movedRight += 1
    # LEFT
    elif random > stay:
        xPos[cell][time+1] = xPos[cell][time] - 1
        print("moved left")
        movedLeft += 1
    # update occupancy based on cell movement
    if random > stay:
        occupied[yPos[cell][time]][xPos[cell][time]] -= 1
        occupied[yPos[cell][time+1]][xPos[cell][time+1]] += 1
    if time % 1 == 0:
        print("fib: " + str(fib[yPos[cell][time-1]][xPos[cell][time-1]]))
        print("vegf: " + str(vegf[yPos[cell][time-1]][xPos[cell][time-1]]))
        print("pro: " + str(pro[yPos[cell][time-1]][xPos[cell][time-1]]))
        print("moves up: " + str(movedUp))
        print("moves down: = " + str(movedDown))
        print("moves right: = " + str(movedRight))
        print("moves left: " + str(movedLeft))

    return xPos, yPos, occupied
