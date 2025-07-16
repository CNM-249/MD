# Local Dynamical Activity Map (YZ Plane)

Compute local dynamical activity maps from molecular dynamics (MD) trajectories, focusing on the 2D plane. This helps visualize spatial variations in molecular mobility.

## Description

This script calculates the mean squared displacement (MSD) of selected atoms between consecutive frames in a trajectory, binned onto a 2D grid in the YZ plane. It can be used to detect regions of high or low molecular motion, which is valuable for studying interfaces, confined fluids, or local dynamics in complex systems.

## Features

- Analyzes MD trajectories in YZ plane
- Uses periodic boundary conditions (minimum image convention)
- Outputs a heatmap visualization of local dynamical activity
- Compatible with GROMACS trajectories via MDAnalysis

## Requirements

- Python ≥ 3.7
- MDAnalysis
- numpy
- matplotlib
- tqdm

Install dependencies (for example, via pip):

## Usage

Update these variables in the script:

gro_file = "xxx.gro"    # Replace with your .gro structure file
xtc_file = "xxx.xtc"    # Replace with your trajectory file
start_frame = 1         # Starting frame index
n_frames = 100          # Number of frames to process
n_bins_y = 600          # Number of grid bins in Y direction
n_bins_z = 340          # Number of grid bins in Z direction

python LMSD.py
