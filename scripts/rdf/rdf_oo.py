"""
rdf_oo.py

Calculate the oxygen–oxygen radial distribution function (RDF)
for bulk water molecular dynamics trajectories using MDAnalysis.

Developed as part of a project.
"""

# ==========================================================
# Import libraries
# ==========================================================

import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from tqdm import tqdm

# ==========================================================
# Load trajectory
# ==========================================================

GRO_FILE = "Bulk-Water/300k/trial.gro"
XTC_FILE = "Bulk-Water/300k/trial.xtc"

u = mda.Universe(GRO_FILE, XTC_FILE)

# Select oxygen atoms
oxygen = u.select_atoms("resname SOL and name OW")

# ==========================================================
# Analysis parameters
# ==========================================================

R_MAX = 10.0          # Maximum RDF distance (Å)
N_BINS = 200          # Number of histogram bins
FRAME_STEP = 1000     # Analyze every 1000th frame

dr = R_MAX / N_BINS
bins = np.linspace(0.0, R_MAX, N_BINS + 1)

rdf_hist = np.zeros(N_BINS)

n_atoms = len(oxygen)
n_frames = 0

# ==========================================================
# RDF calculation
# ==========================================================

for ts in tqdm(u.trajectory[::FRAME_STEP], desc="Calculating O-O RDF"):

    n_frames += 1

    positions = oxygen.positions
    box = u.dimensions

    # Compute all pairwise distances with periodic boundary conditions
    distances = distance_array(positions, positions, box=box)

    # Keep only unique atom pairs (upper triangular matrix)
    upper = np.triu_indices(n_atoms, k=1)
    distances = distances[upper]

    # Accumulate histogram
    hist, _ = np.histogram(distances, bins=bins)
    rdf_hist += hist

# ==========================================================
# Normalization
# ==========================================================

# Bin centers
r = 0.5 * (bins[:-1] + bins[1:])

# Volume of spherical shells
shell_volume = 4.0 * np.pi * r**2 * dr

# Average simulation box volume
volumes = []

for ts in u.trajectory[::FRAME_STEP]:
    lx, ly, lz = u.dimensions[:3]
    volumes.append(lx * ly * lz)

average_volume = np.mean(volumes)

# Number density (Å⁻³)
rho = n_atoms / average_volume

# Normalize RDF
normalization = n_frames * n_atoms * rho

g_r = rdf_hist / (shell_volume * normalization)

# Correct for counting only unique atom pairs
g_r *= 2.0

# ==========================================================
# Save results
# ==========================================================

rdf_df = pd.DataFrame({
    "r (Å)": r,
    "g(r)": g_r
})
rdf_df.to_csv("OUTPUT_FILE", index=False)
