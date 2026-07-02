"""
msd_parallel.py

Compute the mean squared displacement (MSD) of water oxygen atoms
from molecular dynamics trajectories using multiple time-origin
averaging and parallel processing.

The script averages the MSD over multiple independent time origins
to improve statistical sampling. The input trajectory must be
unwrapped so that molecular displacements are continuous across
periodic boundaries.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis import Universe
from joblib import Parallel, delayed
from tqdm import tqdm

GRO_FILE = "md.gro"
XTC_FILE = "md.xtc"
OUTPUT_FILE = "msd.csv"

MSD_TIME = 50.0       # ps
FRAME_INTERVAL = 1000
N_CORES = -1

u = mda.Universe(GRO_FILE, XTC_FILE)
dt = u.trajectory.dt
n_frames = len(u.trajectory)

frame_lag = int(MSD_TIME / dt)
if frame_lag >= n_frames:
    raise ValueError("MSD_TIME is longer than the trajectory duration.")
frame_range = [(start, start + frame_lag)
    for start in range(0, n_frames - frame_lag, FRAME_INTERVAL)]

def calculate_msd(gro_file, xtc_file, start_frame, end_frame):
    u = Universe(gro_file, xtc_file)
    oxygen = u.select_atoms("name OW")
    u.trajectory[start_frame]
    initial_positions = oxygen.positions.copy()
    
    msd = []
    for ts in u.trajectory[start_frame:end_frame]:
        displacement = oxygen.positions - initial_positions
        squared_displacement = np.sum(displacement**2,axis=1)
        msd.append(np.mean(squared_displacement))
    return msd

msd_results = Parallel(n_jobs=N_CORES)(delayed(calculate_msd)(GRO_FILE,XTC_FILE,start,stop) for start, stop in tqdm(frame_range,desc="MSD windows"))

df = pd.DataFrame(msd_results)
mean_msd = df.mean(axis=0)
std_msd = df.std(axis=0)
time = np.arange(frame_lag) * dt
results = pd.DataFrame({
    "Time (ps)": time,
    "MSD (Å²)": mean_msd,
    "STD (Å²)": std_msd})

results.to_csv(OUTPUT_FILE, index=False)
