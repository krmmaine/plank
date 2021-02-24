def move(cell, time, stay, left, right, up, random, yPos, xPos, occupied, fib, vegf, pro, movement, file, T):

    # DOWN
    if random > stay + left + right + up:
        yPos[cell][time+1] = yPos[cell][time] + 1           # don't need to update xpos because I assumed in the beginning that xpos[cell][time+1] = xpos[cell][time]
        if yPos[cell][time + 1] > 0:                        # need time + 1 because the cell just moved down so it is no longer in the capillary
            movement.append(0)
            if len(movement) > 1:
                if movement[len(movement)-2] == 1:
                    print("backtrack")
                    file.write("\n\nBACKTRACK: time = " + str(time) + "\n")
                    file.write("moved down\n\n")

                    file.write("substrate concentrations to the right: \n")
                    file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                    file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                    file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] + 1]) + "\n\n")

                    file.write("substrate concentrations to the left: \n")
                    file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                    file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                    file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] - 1]) + "\n\n")

                    file.write("substrate concentrations to the up: \n")
                    file.write("fib: " + str(fib[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                    file.write("vegf: " + str(vegf[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                    file.write("pro: " + str(pro[yPos[cell][time] - 1][xPos[cell][time]]) + "\n\n")

                    file.write("substrate concentrations to the down: \n")
                    file.write("fib: " + str(fib[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                    file.write("vegf: " + str(vegf[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                    file.write("pro: " + str(pro[yPos[cell][time] + 1][xPos[cell][time]]) + "\n\n")

                    file.write("prob down = " + str(1 - left - right - up - stay) + "\n")
                    file.write("prob up = " + str(up) + "\n")
                    file.write("prob right = " + str(right) + "\n")
                    file.write("prob left = " + str(left) + "\n")
                    file.write("prob stay = " + str(stay) + "\n")

                    file.write("left chemoattractance" + str(T[0]) + "\n")
                    file.write("right chemoattractance" + str(T[1]) + "\n")
                    file.write("up chemoattractance" + str(T[2]) + "\n")
                    file.write("down chemoattractance" + str(T[3]) + "\n")

    # UP
    elif random > stay + left + right:
        yPos[cell][time+1] = yPos[cell][time] - 1
        if yPos[cell][time] > 0:
            movement.append(1)
            if movement[len(movement) - 2] == 0:
                print("backtrack")
                file.write("\n\nBACKTRACK: time = " + str(time) + "\n")
                file.write("moved up\n\n")

                file.write("substrate concentrations to the right: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] + 1]) + "\n\n")

                file.write("substrate concentrations to the left: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] - 1]) + "\n\n")

                file.write("substrate concentrations to the up: \n")
                file.write("fib: " + str(fib[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] - 1][xPos[cell][time]]) + "\n\n")

                file.write("substrate concentrations to the down: \n")
                file.write("fib: " + str(fib[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] + 1][xPos[cell][time]]) + "\n\n")

                file.write("prob down = " + str(1 - left - right - up - stay) + "\n")
                file.write("prob up = " + str(up) + "\n")
                file.write("prob right = " + str(right) + "\n")
                file.write("prob left = " + str(left) + "\n")
                file.write("prob stay = " + str(stay) + "\n")

                file.write("left chemoattractance" + str(T[0]) + "\n")
                file.write("right chemoattractance" + str(T[1]) + "\n")
                file.write("up chemoattractance" + str(T[2]) + "\n")
                file.write("down chemoattractance" + str(T[3]) + "\n")

    # RIGHT
    elif random > stay + left:
        xPos[cell][time+1] = xPos[cell][time] + 1
        if yPos[cell][time] > 0:
            movement.append(2)
            if movement[len(movement) - 2] == 3:
                print("backtrack")
                file.write("\n\nBACKTRACK: time = " + str(time) + "\n")
                file.write("moved right\n\n")

                file.write("substrate concentrations to the right: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] + 1]) + "\n\n")

                file.write("substrate concentrations to the left: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] - 1]) + "\n\n")

                file.write("substrate concentrations to the up: \n")
                file.write("fib: " + str(fib[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] - 1][xPos[cell][time]]) + "\n\n")

                file.write("substrate concentrations to the down: \n")
                file.write("fib: " + str(fib[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] + 1][xPos[cell][time]]) + "\n\n")

                file.write("prob down = " + str(1 - left - right - up - stay) + "\n")
                file.write("prob up = " + str(up) + "\n")
                file.write("prob right = " + str(right) + "\n")
                file.write("prob left = " + str(left) + "\n")
                file.write("prob stay = " + str(stay) + "\n")

                file.write("left chemoattractance" + str(T[0]) + "\n")
                file.write("right chemoattractance" + str(T[1]) + "\n")
                file.write("up chemoattractance" + str(T[2]) + "\n")
                file.write("down chemoattractance" + str(T[3]) + "\n")

    # LEFT
    elif random > stay:
        xPos[cell][time+1] = xPos[cell][time] - 1
        if yPos[cell][time] > 0:
            movement.append(3)
            if movement[len(movement) - 2] == 2:
                print("backtrack")
                file.write("\n\nBACKTRACK: time = " + str(time) + "\n")
                file.write("moved left\n\n")

                file.write("substrate concentrations to the right: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] + 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] + 1]) + "\n\n")

                file.write("substrate concentrations to the left: \n")
                file.write("fib: " + str(fib[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time]][xPos[cell][time] - 1]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time]][xPos[cell][time] - 1]) + "\n\n")

                file.write("substrate concentrations to the up: \n")
                file.write("fib: " + str(fib[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] - 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] - 1][xPos[cell][time]]) + "\n\n")

                file.write("substrate concentrations to the down: \n")
                file.write("fib: " + str(fib[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("vegf: " + str(vegf[yPos[cell][time] + 1][xPos[cell][time]]) + "\n")
                file.write("pro: " + str(pro[yPos[cell][time] + 1][xPos[cell][time]]) + "\n\n")

                file.write("prob down = " + str(1 - left - right - up - stay) + "\n")
                file.write("prob up = " + str(up) + "\n")
                file.write("prob right = " + str(right) + "\n")
                file.write("prob left = " + str(left) + "\n")
                file.write("prob stay = " + str(stay) + "\n")

                file.write("left chemoattractance" + str(T[0]) + "\n")
                file.write("right chemoattractance" + str(T[1]) + "\n")
                file.write("up chemoattractance" + str(T[2]) + "\n")
                file.write("down chemoattractance" + str(T[3]) + "\n")

    # update occupancy based on cell movement
    if random > stay:
        occupied[yPos[cell][time]][xPos[cell][time]] -= 1
        occupied[yPos[cell][time+1]][xPos[cell][time+1]] += 1

    return file
