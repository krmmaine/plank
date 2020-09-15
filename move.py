def move(cell, time, stay, left, right, down, random, yPos, xPos, occupied, fib, vegf, pro, movement):

    file = open("tracking_backtracks.txt", "w")

    # UP
    if random > stay + left + right + down:
        yPos[cell][time+1] = yPos[cell][time] + 1
        if yPos[cell][time] > 0:
            movement.append(0)
            if len(movement) > 1:
                if movement[len(movement)-2] == 1:
                    print("backtrack")
                    file.write("BACKTRACK: time = " + str(time))
                    file.write("moved up")
                    file.write("fib: " + str(fib[yPos[cell][time+1]][xPos[cell][time]]))
                    file.write("vegf: " + str(vegf[yPos[cell][time+1]][xPos[cell][time]]))
                    file.write("pro: " + str(pro[yPos[cell][time+1]][xPos[cell][time]]))
                    file.write("prob up = " + str(1-left-right-down-stay))
                    file.write("prob down = " + str(down))
                    file.write("prob right = " + str(right))
                    file.write("prob left = " + str(left))

    # DOWN
    elif random > stay + left + right:
        yPos[cell][time+1] = yPos[cell][time] - 1
        if yPos[cell][time] > 0:
            movement.append(1)
            if movement[len(movement) - 2] == 0:
                print("backtrack")
                file.write("BACKTRACK: time = " + str(time))
                file.write("moved down")
                file.write("fib: " + str(fib[yPos[cell][time + 1]][xPos[cell][time]]))
                file.write("vegf: " + str(vegf[yPos[cell][time + 1]][xPos[cell][time]]))
                file.write("pro: " + str(pro[yPos[cell][time + 1]][xPos[cell][time]]))
                file.write("prob up = " + str(1 - left - right - down - stay))
                file.write("prob down = " + str(down))
                file.write("prob right = " + str(right))
                file.write("prob left = " + str(left))

    # RIGHT
    elif random > stay + left:
        xPos[cell][time+1] = xPos[cell][time] + 1
        if yPos[cell][time] > 0:
            movement.append(2)
            if movement[len(movement) - 2] == 3:
                print("backtrack")
                file.write("BACKTRACK: time = " + str(time))
                file.write("moved right")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("prob up = " + str(1 - left - right - down - stay))
                file.write("prob down = " + str(down))
                file.write("prob right = " + str(right))
                file.write("prob left = " + str(left))

    # LEFT
    elif random > stay:
        xPos[cell][time+1] = xPos[cell][time] - 1
        if yPos[cell][time] > 0:
            movement.append(3)
            if movement[len(movement) - 2] == 2:
                print("backtrack")
                file.write("BACKTRACK: time = " + str(time))
                file.write("moved left")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time + 1]]))
                file.write("prob up = " + str(1 - left - right - down - stay))
                file.write("prob down = " + str(down))
                file.write("prob right = " + str(right))
                file.write("prob left = " + str(left))

    # update occupancy based on cell movement
    if random > stay:
        occupied[yPos[cell][time]][xPos[cell][time]] -= 1
        occupied[yPos[cell][time+1]][xPos[cell][time+1]] += 1

    file.close()

    return movement
