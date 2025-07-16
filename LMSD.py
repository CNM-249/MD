import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Set input files and parameters
gro_file = "xxx.gro"
xtc_file = "xxx.xtc"
start_frame = 1
n_frames = 100  # Total frames to analyze

n_bins_y = 600
n_bins_z = 340

# Load trajectory
u = mda.Universe(gro_file, xtc_file)
ow_atoms = u.select_atoms("resname CO2 and name C")

# Get box size (assuming fixed box dimensions)
u.trajectory[start_frame]
box = u.trajectory.ts.dimensions[:3]
box_y, box_z = box[1], box[2]

# Define grid edges
y_edges = np.linspace(0, box_y, n_bins_y + 1)
z_edges = np.linspace(0, box_z, n_bins_z + 1)

# Initialize cumulative squared displacements and counters
msd_map = np.zeros((n_bins_y, n_bins_z))
count_map = np.zeros((n_bins_y, n_bins_z))

# Main loop over frames
for i in tqdm(range(start_frame, start_frame + n_frames - 1)):
    u.trajectory[i]
    pos_i = ow_atoms.positions.copy()
    box_i = u.trajectory.ts.dimensions[:3]

    u.trajectory[i + 1]
    pos_ip1 = ow_atoms.positions.copy()
    box_ip1 = u.trajectory.ts.dimensions[:3]

    # Minimum image convention correction (only in y and z)
    delta_yz = pos_ip1[:, 1:3] - pos_i[:, 1:3]
    delta_yz -= np.round(delta_yz / box_i[1:3]) * box_i[1:3]
    dr2 = np.sum(delta_yz**2, axis=1)

    # Determine which grid cell each atom belongs to (using frame i)
    y_idx = np.digitize(pos_i[:, 1], y_edges) - 1
    z_idx = np.digitize(pos_i[:, 2], z_edges) - 1

    # Exclude out-of-bounds indices
    valid = (y_idx >= 0) & (y_idx < n_bins_y) & (z_idx >= 0) & (z_idx < n_bins_z)
    y_idx = y_idx[valid]
    z_idx = z_idx[valid]
    dr2 = dr2[valid]

    for y, z, d2 in zip(y_idx, z_idx, dr2):
        msd_map[y, z] += d2
        count_map[y, z] += 1

# Compute average squared displacement
with np.errstate(divide='ignore', invalid='ignore'):
    avg_msd_map = np.where(count_map > 0, msd_map / count_map, 0)

# Plotting
plt.figure(figsize=(7.5, 6))

# Note: extent = [left, right, bottom, top], corresponding to y and z ranges
extent = [0, box_y, 0, box_z]  # x-axis = y, y-axis = z

# Transpose avg_msd_map so rows correspond to z, columns to y
plt.imshow(avg_msd_map.T, origin='lower', extent=extent, aspect='auto',
           cmap='viridis', vmin=0, vmax=8)  # adjust vmin/vmax if needed

# Add labels and colorbar
plt.xlabel("Y (nm)")
plt.ylabel("Z (nm)")
plt.title("Local Dynamical Activity Map (YZ Plane)")
cbar = plt.colorbar()
cbar.set_label("Average $\Delta r_{yz}^2$ (nm$^2$)")
plt.tight_layout()
plt.savefig("yz_LMSD_map.png", dpi=300)
plt.show()
