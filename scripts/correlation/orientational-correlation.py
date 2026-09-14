"""
O-H rotational correlation function of water.

This script calculates the first-, second-, and third-rank rotational
correlation functions (P1, P2, and P3) of O-H bond vectors.

The same correlation-function calculation can also be applied to the
molecular dipole vector. To calculate the dipole rotational correlation,
replace the O-H vector construction with the normalized molecular
dipole vectors while keeping the correlation calculation unchanged.
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
import scipy.optimize
import warnings
from scipy.spatial.distance import cdist, distance_matrix, distance
from MDAnalysis.analysis import distances
from MDAnalysis.lib._cutil import _in2d
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_bonds
from MDAnalysis.core.groups import requires
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
from tqdm import tqdm
from joblib import Parallel, delayed
import seaborn as sns
np.set_printoptions(threshold=np.inf)
np.set_printoptions(threshold=10.80)
font = {
    'family': 'serif',
    'color': 'darkred',
    'weight': 'normal',
    'size': 16
}


def dipole(pos, box):
    oxygen = pos[::3]
    hydrogen1 = pos[1::3]
    hydrogen2 = pos[2::3]
    vec_1 = oxygen - hydrogen1
    vec_2 = oxygen - hydrogen2
    vec_1 -= np.rint(vec_1 / box) * box
    vec_2 -= np.rint(vec_2 / box) * box
    dipole_vector = vec_1 + vec_2
    norm_dipole_vector = (dipole_vector / np.linalg.norm(dipole_vector, axis=1, keepdims=True))
    return norm_dipole_vector
'''
def c2_single_particle(arr1, arr2):
    normalized_arr1 = (arr1 / np.linalg.norm(arr1, axis=1, keepdims=True))
    normalized_arr2 = (arr2 / np.linalg.norm(arr2, axis=1, keepdims=True))
    dot_products = np.dot(normalized_arr1,normalized_arr2.T).diagonal()
    p2 = 0.5 * ((3 * (dot_products ** 2)) - 1)
    return p2

def oh_compute_correlation_singleparticle(frame_pair):
    u.trajectory[frame_pair[0]]
    ini_oh_vector = oh_vec(water_selection.positions,box=u.dimensions[0:3])
    u.trajectory[frame_pair[1]]
    final_oh_vector = oh_vec(water_selection.positions,box=u.dimensions[0:3])
    corr = c2_single_particle(ini_oh_vector,final_oh_vector).flatten()
    return sum(corr), len(corr)
'''
def oh_vec(pos, box=None):
    oxygen = pos[::3]
    hydrogen1 = pos[1::3]
    hydrogen2 = pos[2::3]
    delta_OH1 = oxygen - hydrogen1
    delta_OH2 = oxygen - hydrogen2
    if box is not None:
        delta_OH1 -= np.rint(delta_OH1 / box) * box
        delta_OH2 -= np.rint(delta_OH2 / box) * box
    norm_OH1 = np.linalg.norm(delta_OH1,axis=1,keepdims=True)
    norm_OH2 = np.linalg.norm(delta_OH2,axis=1,keepdims=True)
    normalised_vector_1 = delta_OH1 / norm_OH1
    normalised_vector_2 = delta_OH2 / norm_OH2
    combined_vectors = np.concatenate((normalised_vector_1, normalised_vector_2),axis=0)
    return combined_vectors

def normalize_vector(arr):
    return arr / np.linalg.norm(arr)

def compute_oh_vectors(pos, box):
    pos = pos.reshape(n_waters, 3, 3)
    O = pos[:, 0]
    H1 = pos[:, 1]
    H2 = pos[:, 2]
    v1 = O - H1
    v2 = O - H2
    v1 -= np.rint(v1 / box) * box
    v2 -= np.rint(v2 / box) * box
    v1 /= np.linalg.norm(v1, axis=1, keepdims=True)
    v2 /= np.linalg.norm(v2, axis=1, keepdims=True)
    return np.concatenate((v1, v2), axis=0 )

gro = 'dynamic.gro'
xtc = 'dynamic.xtc'
correlation_time = 15
u = mda.Universe(gro, xtc)
dt = u.trajectory.dt
ow = u.select_atoms('name OW')
resid_confined = ow.resids
water_selection = u.select_atoms(f"resname SOL and resid {' '.join(map(str, resid_confined))}")
n_atoms = len(water_selection)
n_waters = n_atoms // 3
total_frames = int(len(u.trajectory) * 0.9)
frame_lag = int(correlation_time / dt)
all_vectors = np.zeros((total_frames, 2 * n_waters, 3),dtype=np.float32)
start = time.time()

for i, ts in enumerate(tqdm(u.trajectory[0:total_frames],desc="Calculating O-H vectors")):
    box = ts.dimensions[:3]
    pos = water_selection.positions.copy()
    all_vectors[i] = compute_oh_vectors(pos,box)
end = time.time()
print(f"Vector loading time: {end - start:.2f} s")

def compute_lag(lag):
    valid = total_frames - lag
    v0 = all_vectors[:valid]
    vt = all_vectors[lag:lag + valid]
    dot = np.sum(v0 * vt, axis=2)
    p1 = dot
    p2 = 0.5 * (3 * dot**2 - 1)
    p3 = 0.5 * (5 * dot**3 - 3 * dot)
    return (p1.mean(),p2.mean(),p3.mean() )

corr = Parallel(n_jobs=55)(delayed(compute_lag)(lag)for lag in tqdm(range(frame_lag),desc="Calculating correlation"))
corr = np.array(corr)
corr /= corr[0]
time_axis = np.arange(frame_lag) * dt
df_single = pd.DataFrame(corr)
df_single.columns = ['p1','p2','p3']
df_single['time'] = time_axis
df_single = df_single[['time', 'p1', 'p2', 'p3']]
df_single.to_csv('tcf-15ps-p1-p2-p3.csv',index=False)
