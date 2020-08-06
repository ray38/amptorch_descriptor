import numpy as np
from amptorch_descriptor.BP_symmetry_function import BPSymmetryFunction
from amptorch_descriptor.dataset import AMPTorchDataset
from ase.io.trajectory import Trajectory
from ase.io import read

elements = ["H","O"]
Gs = {"G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0] * 4},\
      "G4": {"etas": [0.005], "zetas": [1.0], "gammas": [-1.0, 1.0]},\
      "cutoff": 6.5}

small = read('./small/water.extxyz', index=':')
trajectories = [small]

descriptor = BPSymmetryFunction(Gs = Gs, elements = elements)

training_data = AMPTorchDataset(trajectories, descriptor)