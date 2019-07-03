# Set of python codes to compute covariance of COM mapped protein trajectory


Requirements:
MDAnalysis
numpy

## traj_routines.py

	routines used by other python scripts to align trajectories

## cgcom_align_traj.py

	python script to read set of atomistic trajectories and create a COM mapped and aligned trajectory

## covar_10ns_windows.py

	python script to read coarse-grained trajectory and create covariance matrices and average structure for chunks of trajectory

## avg_covar.py

	python script to read covariance matrices and average structures from previous output, align average structures, and rotate/average covariance matrices
