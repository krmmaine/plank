def move(cell, time, stay, left, right, down, random, yPos, xPos, occupied):

    # UP
    if random > stay + left + right + down:
        yPos[cell][time+1] = yPos[cell][time] + 1
    # DOWN
    elif random > stay + left + right:
        yPos[cell][time+1] = yPos[cell][time] - 1
    # RIGHT
    elif random > stay + left:
        xPos[cell][time+1] = xPos[cell][time] + 1
    # LEFT
    elif random > stay:
        xPos[cell][time+1] = xPos[cell][time] - 1
    # update occupancy based on cell movement
    if random > stay:
        occupied[yPos[cell][time]][xPos[cell][time]] -= 1
        occupied[yPos[cell][time+1]][xPos[cell][time+1]] += 1

    return xPos, yPos, occupied
