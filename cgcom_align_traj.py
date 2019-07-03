import numpy as np
import MDAnalysis as md
import traj_routines as traj

selection = "resid 1:454"

# setup MDAnalysis universe and selection
fileList = []
for j in range(4):
    fileBase = "/Users/martinmccullagh/Documents/research/2019/allostery/IGPS_traj/hk181a/system" + str(j+1) + "/hK181A_mutant.prod."
    for i in range(26,51):
        fileName = fileBase + str(i) + ".dcd"
        fileList.append(fileName)
coord = md.Universe("/Users/martinmccullagh/Documents/research/2019/allostery/IGPS_traj/hk181a/truncated.hK181A_mutant.prmtop",fileList)
print("Number of frames:",coord.trajectory.n_frames)
sel = coord.select_atoms(selection)
# iteratively align and average trajectory
avgPos, alignedPos  = traj.iterative_align_average_com(coord,sel,thresh=1E-10)[:2]
# save average structure
np.savetxt("IGPS_apo_hK181A_1us_com_average_structure.dat",avgPos)

# open trajectory file for writing
ca_sel = coord.select_atoms("name CA")
positions_xyz = md.Writer("IGPS_apo_hK181A_com_positions.dcd",ca_sel.n_atoms)
ca_sel.positions = avgPos
# save a pdb of average structure
ca_sel.write("IGPS_apo_hK181A_1us_com_average_structure.pdb")

# loop through trajectory and compute covariance 
for ts in coord.trajectory:
    ca_sel.positions = alignedPos[ts.frame-1,:,:]
    positions_xyz.write(ca_sel)
positions_xyz.close()

