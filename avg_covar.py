
import numpy as np
import traj_routines as traj

#############################################################################################################
###################################    SUBROUTINES     ######################################################
#############################################################################################################

def rot_covar(covar,rot):
    N3 = covar.shape[0]
    N = N3//3
    for i in range(N):
        for j in range(N):
            covar[i*3:(i+1)*3,j*3:(j+1)*3] = np.dot(covar[i*3:(i+1)*3,j*3:(j+1)*3],rot)
    return covar

#############################################################################################################
###################################    MAIN PROGRAM    ######################################################
#############################################################################################################

# number of windows to read/average
nWindows = 100
# input file name variables
root = "IGPS_holo_window_"
avgEnd = "_average.dat"
covarEnd = "_covar.dat"

# read averages
fileName = root + str(1) + avgEnd
avg1 = np.loadtxt(fileName)
avgTraj = np.empty((nWindows,avg1.shape[0],avg1.shape[1]),dtype=np.float32)
avgTraj[0,:,:] = avg1

for window in range(1,nWindows):

    fileName = root + str(window+1) + avgEnd
    avgTraj[window,:,:] = np.loadtxt(fileName)

# align average structures and keep rotation matrices
avg, avgTraj, rot = traj.iterative_align_coord_array(avgTraj,maxSteps=100,thresh=1E-6)
# write average of averages to file
np.savetxt("IGPS_holo_10ns_window_average_structure.dat",avg)

# read covariances and rotate them
fileName = root + str(1) + covarEnd
covar = rot_covar(np.loadtxt(fileName),rot[0,:,:])
for window in range(1,nWindows):
    fileName = root + str(window+1) + covarEnd
    covar += rot_covar(np.loadtxt(fileName),rot[window,:,:])
covar /= nWindows
# write new covariance
np.savetxt("IGPS_holo_10ns_window_average_covar.dat",covar)

