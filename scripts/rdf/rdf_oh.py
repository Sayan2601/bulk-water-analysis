"""
rdf_oh.py

Compute the intermolecular oxygen–hydrogen (O–H) or hydrogen–oxygen
(H–O) radial distribution function from bulk water molecular dynamics
trajectories using MDAnalysis.

To calculate the H–O RDF, simply swap the atom selections below.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm

# Input files
GRO_FILE = "md.gro"
XTC_FILE = "md.xtc"
OUTPUT_FILE = "rdf_oh.csv"

# Atom selections
REFERENCE_SELECTION = "resname SOL and name OW"
NEIGHBOR_SELECTION = "resname SOL and name HW1 HW2"

# For H(reference)–O(neighbor), replace with:
#
# REFERENCE_SELECTION = "resname SOL and name HW1 HW2"
# NEIGHBOR_SELECTION  = "resname SOL and name OW"
# OUTPUT_FILE = "rdf_ho.csv"

# Analysis parameters
R_MAX = 10.0          # Å
N_BINS = 200
FRAME_STEP = 1

u = mda.Universe(GRO_FILE, XTC_FILE)
reference = u.select_atoms(REFERENCE_SELECTION)
neighbor = u.select_atoms(NEIGHBOR_SELECTION)
n_reference = len(reference)
n_neighbor = len(neighbor)
bins = np.linspace(0.0, R_MAX, N_BINS + 1)
dr = bins[1] - bins[0]
rdf_hist = np.zeros(N_BINS)
n_frames = 0
volumes = []

for ts in tqdm(u.trajectory[::FRAME_STEP], desc="Calculating RDF"):
    n_frames += 1
    distances = distance_array(reference.positions,neighbor.positions,box=u.dimensions)
    same_residue = (reference.resids[:, None] ==neighbor.resids[None, :])
    distances = distances[~same_residue]
    hist, _ = np.histogram(distances, bins=bins)
    rdf_hist += hist
    volumes.append(np.prod(u.dimensions[:3]))

average_volume = np.mean(volumes)
rho_neighbor = n_neighbor / average_volume
r = 0.5 * (bins[:-1] + bins[1:])
shell_volume = 4.0 * np.pi * r**2 * dr
g_r = rdf_hist / (n_frames * n_reference *rho_neighbor *shell_volume)
rdf_df = pd.DataFrame({"r (Å)": r, "g(r)": g_r})
rdf_df.to_csv(OUTPUT_FILE, index=False)
