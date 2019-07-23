# load libraries
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

kb = 1.9872036E-3 # kcal/mol/K

def iterative_align_average(coord,selGroup,frameStart=0,frameStop=-1,deltaFrame=1,maxSteps=25,thresh=0.001):
    if frameStop < 0:
        frameStop = coord.trajectory.n_frames + frameStop + 1
    nFrames= int((frameStop-frameStart)//deltaFrame)
    # create numpy array of aligned positions
    alignedPos = np.empty((nFrames,selGroup.n_atoms,3),dtype=np.float64)
    # create numpy array of rotation matrices
    rotMat = np.empty((nFrames,3,3),dtype=np.float64)
    #generate an initial average by aligning to first frame
    avg = np.zeros((selGroup.n_atoms,3),dtype=np.float64)
    frameCount = 0
    for ts in coord.trajectory[frameStart:frameStop:deltaFrame]:
        selGroup.translate(-selGroup.center_of_mass())
        if frameCount == 0:
            ref = np.copy(selGroup.positions)
        else:
            R = align.rotation_matrix(selGroup.positions, ref)[0]
            rotMat[frameCount,:,:] = np.copy(R)
            selGroup.rotate(R)
        avg += selGroup.positions
        alignedPos[frameCount,:,:] = selGroup.positions
        frameCount += 1
    # finish average
    avg /= nFrames
    # perform iterative alignment and average to converge average
    newAvg = np.zeros((selGroup.n_atoms,3),dtype=np.float64)
    avgRmsd = 2*thresh
    step = 0
    while step<maxSteps and avgRmsd > thresh:
        newAvg = 0.0
        frameCount = 0
        for ts in coord.trajectory[frameStart:frameStop:deltaFrame]:
            alignedPos[frameCount,:,:] -= np.mean(alignedPos[frameCount,:,:],axis=0)
            R = align.rotation_matrix(alignedPos[frameCount,:,:], avg)[0]
            rotMat[frameCount,:,:]  = np.dot(rotMat[frameCount,:,:],R)
            alignedPos[frameCount,:,:] = np.dot(alignedPos[frameCount,:,:],R.T)
            newAvg += alignedPos[frameCount,:,:]
            frameCount += 1
        # finish average
        newAvg /= nFrames
        avgRmsd = rmsd(avg,newAvg,center=False,superposition=False)
        avg = np.copy(newAvg)
        step += 1
        print(step, avgRmsd)
    return avg, alignedPos, rotMat


def iterative_align_average_com(coord,selGroup,frameStart=0,frameStop=-1,deltaFrame=1,maxSteps=25,thresh=0.001):
    if frameStop < 0:
        frameStop = coord.trajectory.n_frames + frameStop + 1
    nFrames= int((frameStop-frameStart)//deltaFrame)
    # create numpy array of aligned positions
    alignedPos = np.empty((nFrames,selGroup.n_residues,3),dtype=np.float64)
    #generate an initial average by aligning to first frame
    avg = np.zeros((selGroup.n_residues,3),dtype=np.float64)
    comPos = np.empty((selGroup.n_residues,3),dtype=np.float64)
    frameCount = 0
    for ts in coord.trajectory[frameStart:frameStop:deltaFrame]:
        selGroup.translate(-selGroup.center_of_mass())
        for i, resid in enumerate(np.unique(selGroup.resids)):
            residSel = "resid " + str(resid)
            comPos[i,:] = selGroup.select_atoms(residSel).center_of_mass()            
        if frameCount == 0:
            ref = np.copy(comPos)
        else:
            R = align.rotation_matrix(comPos, ref)[0]
            comPos = np.dot(comPos,R.T)
        avg += comPos
        alignedPos[frameCount,:,:] = comPos
        frameCount += 1
    # finish average
    avg /= nFrames
    # perform iterative alignment and average to converge average
    newAvg = np.zeros((selGroup.n_residues,3),dtype=np.float64)
    avgRmsd = 2*thresh
    step = 0
    while step<maxSteps and avgRmsd > thresh:
        newAvg = 0.0
        frameCount = 0
        for ts in coord.trajectory[frameStart:frameStop:deltaFrame]:
            alignedPos[frameCount,:,:] -= np.mean(alignedPos[frameCount,:,:],axis=0)
            R = align.rotation_matrix(alignedPos[frameCount,:,:], avg)[0]
            alignedPos[frameCount,:,:] = np.dot(alignedPos[frameCount,:,:],R.T)
            newAvg += alignedPos[frameCount,:,:]
            frameCount += 1
        # finish average
        newAvg /= nFrames
        avgRmsd = rmsd(avg,newAvg,center=False,superposition=False)
        avg = np.copy(newAvg)
        step += 1
        print(step, avgRmsd)
    return avg, alignedPos

def iterative_align_coord_array(coord,maxSteps=20,thresh=0.001):
    nFrames= coord.shape[0]
    nAtoms = coord.shape[1]
    # save rotation matrices
    rot = np.empty((nFrames,3,3),dtype=float)
    #generate an initial average by aligning to first frame
    avg = np.zeros((nAtoms,3),dtype=np.float64)
    for ts in range(nFrames):
        coord[ts,:,:] -= np.mean(coord[ts,:,:],axis=0)
        if ts == 0:
            ref = np.copy(coord[ts,:,:])
            rot[ts,:,:] = np.identity(3)
        else:
            R = align.rotation_matrix(coord[ts,:,:], ref)[0]
            rot[ts,:,:] = R.T
            coord[ts,:,:] = np.dot(coord[ts,:,:],R.T)
        avg += coord[ts,:,:]
    # finish average
    avg /= nFrames
    # perform iterative alignment and average to converge average
    newAvg = np.zeros((nAtoms,3),dtype=np.float64)
    avgRmsd = 2*thresh
    step = 0
    while step<maxSteps and avgRmsd > thresh:
        newAvg = 0.0
        for ts in range(nFrames):
            coord[ts,:,:] -= np.mean(coord[ts,:,:],axis=0)
            R = align.rotation_matrix(coord[ts,:,:], avg)[0]
            rot[ts,:,:] = np.dot(rot[ts,:,:],R.T)
            coord[ts,:,:] = np.dot(coord[ts,:,:],R.T)
            newAvg += coord[ts,:,:]
        # finish average
        newAvg /= nFrames
        avgRmsd = rmsd(avg,newAvg,center=False,superposition=False)
        avg = np.copy(newAvg)
        step += 1
        print(step, avgRmsd)
    return avg, coord, rot

