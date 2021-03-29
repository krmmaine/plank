import numpy as np
import csv

# File takes in the cell positions at each time (x and y)
# Calculates the FMIx, FMIy, Directness, Euclidean distance, Velocity

def data_output(celltrackingvector,xlength,xsteps,totalTime,numTimesteps):
    with open('migration_data.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['Cell Number', 'Accumulated Distance', 'Euclidian Distance', 'Directness', 'FMI_y', 'FMI_x',
                             'Accumulated Distance (mm)', 'Euclidian Distance (mm)', 'Euclidian Velocity (mm/hr)',
                             'Accumulated Velocity (mm/hr) ', 'Time (hr)'])

        for cell in range(len(celltrackingvector)):
            # Don't measure empty vectors
            if len(celltrackingvector[cell][0]) == 0:
                continue

            # Loop through data and remove the data when the cells are still in the parent vessel
            for i in range(len(celltrackingvector[cell][0])):
                if celltrackingvector[cell][2][i] != 0:
                    index = i
                    break
            del celltrackingvector[cell][0][:index]
            del celltrackingvector[cell][1][:index]
            del celltrackingvector[cell][2][:index]

            distances = [] # Initialize vector to collect distance moved each timestep
            for time in range(1,len(celltrackingvector[cell][0])):
                # Calculate the distance at each time step using the pythagorean formula
                distances.append(np.sqrt((celltrackingvector[cell][1][time] - celltrackingvector[cell][1][time-1])**2 +
                                         (celltrackingvector[cell][2][time] - celltrackingvector[cell][2][time-1])**2))

            # Calculate outputs
            dist_accum = sum(distances)
            dist_euclid = np.sqrt((celltrackingvector[cell][1][-1] - celltrackingvector[cell][1][0])**2 +
                                         (celltrackingvector[cell][2][-1] - celltrackingvector[cell][2][0])**2)
            directness = dist_euclid / dist_accum
            FMI_y = (celltrackingvector[cell][2][-1] - celltrackingvector[cell][2][0]) / dist_accum
            FMI_x = (celltrackingvector[cell][1][-1] - celltrackingvector[cell][1][0]) / dist_accum

            dist_accum_nondimensionalized = dist_accum * (xlength/xsteps)
            dist_euclid_nondimensionalized = dist_euclid * (xlength/xsteps)

            time_nondimensionalized = (celltrackingvector[cell][0][-1] - celltrackingvector[cell][0][0]) * (totalTime/numTimesteps)
            velocity_euclid_nondimensionalized = dist_euclid_nondimensionalized / time_nondimensionalized
            velocity_accum_nondimensionalized = dist_accum_nondimensionalized / time_nondimensionalized

            filewriter.writerow([cell,dist_accum,dist_euclid,directness,FMI_y,FMI_x,
                                 dist_accum_nondimensionalized,dist_euclid_nondimensionalized, velocity_euclid_nondimensionalized,
                                 velocity_accum_nondimensionalized, time_nondimensionalized])

        #print(celltrackingvector)
    return
