import numpy as np
from .BP_symmetry_function import BP_symmetry_function
from .dataset import AMPTorchDataset
from .util import Cosine_cutoff

elements = ["H","C","O"]
Gs = {"G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0] * 4},\
      "G4": {"etas": [0.005], "zetas": [1.0], "gammas": [-1.0, 1.0]},\
      "Cutoff": Cosine_cutoff(6.5)}

trajectories = []

descriptor = BP_symmetry_function(Gs = Gs, elements = elements)

training_data = AMPTorchDataset