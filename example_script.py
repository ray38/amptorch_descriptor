import numpy as np
from amptorch_descriptor.BP_symmetry_function import BPSymmetryFunction
from amptorch_descriptor.Atomistic_MCSH import AtomisticMCSH
# from amptorch_descriptor.dataset import AMPTorchDataset
from amptorch_descriptor.descriptor_calculator import DescriptorCalculator
from ase.io.trajectory import Trajectory
from ase.io import read

large = Trajectory('./large/iron_data.traj')
trajectories = [large]

elements = ["H","O","Fe"]



Gs = {"G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0] * 4},\
      "G4": {"etas": [0.005], "zetas": [1.0], "gammas": [-1.0, 1.0]},\
      "cutoff": 6.5}

MCSHs = {   "MCSHs": {   
                        "0": {"groups": [1], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "1": {"groups": [1], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "2": {"groups": [1,2], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "3": {"groups": [1,2,3], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "4": {"groups": [1,2,3,4], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "5": {"groups": [1,2,3,4,5], "sigmas": [0.1, 0.2, 0.3, 0.4]},
                        "6": {"groups": [1,2,3,4,5,6,7], "sigmas": [0.1, 0.2, 0.3, 0.4]}
                  },
            "atom_gaussians": {
                        "H": "./MCSH_potential/H_pseudodensity_6.g",
                        "O": "./MCSH_potential/O_pseudodensity_6.g",
                        "Fe": "./MCSH_potential/Pt_pseudodensity_8.g"
                  },
            "cutoff": 6.5}

# small = read('./small/water.extxyz', index=':')
# trajectories = [small]


# descriptor = BPSymmetryFunction(Gs = Gs, elements = elements)

descriptor = AtomisticMCSH(MCSHs = MCSHs, elements = elements)
descriptor_calculator = DescriptorCalculator(trajectories, descriptor)

descriptor_calculator.calculate_PCA(n_components = 10)