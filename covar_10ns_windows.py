
import numpy as np
import MDAnalysis as md
import traj_routines as traj

nWindows = 100
# coarse-grained trajectory input
topFile = "CHANGEME.pdb"
trajFile = "CHANGEME.dcd"
coord = md.Universe(topFile,trajFile)
sel = coord.select_atoms("all")

deltaSteps = coord.trajectory.n_frames//nWindows

# file name variables
root = "IGPS_holo_window_"
avgEnd = "_average.dat"
covarEnd = "_covar.dat"

for window in range(nWindows):
    # set some window variables
    frameStart = window * deltaSteps
    frameStop = frameStart + deltaSteps
    # compute windows average position and aligned positions
    avgPos, alignedPos = traj.iterative_align_average(coord,sel,frameStart=frameStart,frameStop=frameStop,thresh=1E-10)[:2]
    fileName = root + str(window+1) + avgEnd
    np.savetxt(fileName,avgPos)

    # loop through trajectory and compute covariance 
    covar = np.zeros((3*sel.n_atoms,3*sel.n_atoms),dtype=np.float64)
    for frame in range(deltaSteps):
        covar += np.outer(alignedPos[frame,:,:].reshape(3*sel.n_atoms),alignedPos[frame,:,:].reshape(3*sel.n_atoms))
    # finish covariance
    covar /= deltaSteps
    covar -= np.outer(avgPos.reshape(3*sel.n_atoms),avgPos.reshape(3*sel.n_atoms))
    fileName = root + str(window+1) + covarEnd
    np.savetxt(fileName,covar)

