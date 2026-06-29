"""
rdf_oh.py

Calculates the intermolecular O-H/H-O radial distribution function.
Intramolecular O-H pairs are excluded.

Examples
--------
Oxygen (reference) - Hydrogen:
    REFERENCE_SELECTION = "resname SOL and name OW"
    NEIGHBOR_SELECTION  = "resname SOL and name HW1 HW2"

Hydrogen (reference) - Oxygen:
    REFERENCE_SELECTION = "resname SOL and name HW1 HW2"
    NEIGHBOR_SELECTION  = "resname SOL and name OW"

Developed as part of a  project.
"""

# ==========================================================
# Import libraries
# ==========================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm

# ==========================================================
# Input files
# ==========================================================

GRO_FILE = "md.gro"
XTC_FILE = "md.xtc"

# ==========================================================
# Atom selections
# ==========================================================

# -------- O(reference) - H(neighbour) --------
REFERENCE_SELECTION = "resname SOL and name OW"
NEIGHBOR_SELECTION = "resname SOL and name HW1 HW2"
OUTPUT_FILE = "rdf_oh.csv"
PLOT_TITLE = "Intermolecular O-H RDF"
YLABEL = r"$g_{OH}(r)$"

# ----------------------------------------------------------
# Uncomment the following block for H(reference)-O(neighbour)
#
# REFERENCE_SELECTION = "resname SOL and name HW1 HW2"
# NEIGHBOR_SELECTION = "resname SOL and name OW"
# OUTPUT_FILE = "rdf_ho.csv"
# PLOT_TITLE = "Intermolecular H-O RDF"
# YLABEL = r"$g_{HO}(r)$"
# ----------------------------------------------------------

# ==========================================================
# Parameters
# ==========================================================

R_MAX = 10.0
N_BINS = 200
FRAME_STEP = 100

bins = np.linspace(0, R_MAX, N_BINS + 1)
dr = bins[1] - bins[0]

rdf_hist = np.zeros(N_BINS)

# ==========================================================
# Load trajectory
# ==========================================================

u = mda.Universe(GRO_FILE, XTC_FILE)

reference = u.select_atoms(REFERENCE_SELECTION)
neighbor = u.select_atoms(NEIGHBOR_SELECTION)

n_reference = len(reference)
n_neighbor = len(neighbor)

n_frames = 0

# ==========================================================
# RDF calculation
# ==========================================================

for ts in tqdm(u.trajectory[::FRAME_STEP], desc="Calculating RDF"):

    n_frames += 1

    dists = distance_array(
        reference.positions,
        neighbor.positions,
        box=u.dimensions
    )

    # Remove intramolecular pairs
    same_residue = (
        reference.resids[:, None] ==
        neighbor.resids[None, :]
    )

    dists = dists[~same_residue]

    hist, _ = np.histogram(dists, bins=bins)
    rdf_hist += hist

# ==========================================================
# Normalization
# ==========================================================

volumes = []

for ts in u.trajectory[::FRAME_STEP]:
    lx, ly, lz = u.dimensions[:3]
    volumes.append(lx * ly * lz)

average_volume = np.mean(volumes)

rho_neighbor = n_neighbor / average_volume

r = 0.5 * (bins[:-1] + bins[1:])
shell_volume = 4.0 * np.pi * r**2 * dr

g_r = rdf_hist / (
    n_frames *
    n_reference *
    rho_neighbor *
    shell_volume
)

# ==========================================================
# Save results
# ==========================================================

rdf_df = pd.DataFrame({
    "r (Å)": r,
    "g(r)": g_r
})

rdf_df.to_csv(OUTPUT_FILE, index=False)
