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
                        "0": {"groups": [1], "sigmas": [0.01]},
                        #"1": {"groups": [1], "sigmas": [0.1, 0.2, 0.3]},
                        #"2": {"groups": [1,2], "sigmas": [0.1, 0.2, 0.3]},
                        #"3": {"groups": [1,2,3], "sigmas": [0.1, 0.2, 0.3]},
                        #"4": {"groups": [1,2,3,4], "sigmas": [0.1, 0.2, 0.3]},
                        #"5": {"groups": [1,2,3,4,5], "sigmas": [0.1, 0.2, 0.3]},
                        #"6": {"groups": [1,2,3,4,5,6,7], "sigmas": [0.1, 0.2, 0.3]},
                        #"7": {"groups": [1,2,3,4,5,6,7,8], "sigmas": [0.1, 0.2, 0.3]},
                        #"8": {"groups": [1,2,3,4,5,6,7,8,9,10], "sigmas": [0.1, 0.2, 0.3]},
                        #"9": {"groups": [1,2,3,4,5,6,7,8,9,10,11,12], "sigmas": [0.1, 0.2, 0.3]}
                  },
            "atom_gaussians": {
                        "H": "./MCSH_potential/H_pseudodensity_6.g",
                        "O": "./MCSH_potential/O_pseudodensity_6.g",
                        "Fe": "./MCSH_potential/Pt_pseudodensity_8.g"
                  },
            "cutoff": 3.5,
            "prime_threshold": 1e-6
}

# small = read('./small/water.extxyz', index=':')
# trajectories = [small]

print(large[0].cell)
descriptor = BPSymmetryFunction(Gs = Gs, elements = elements)

#descriptor = AtomisticMCSH(MCSHs = MCSHs, elements = elements)
descriptor_calculator = DescriptorCalculator(trajectories, descriptor)

temp_desc = descriptor_calculator.get_descriptors()
# temp_desc_prime = descriptor_calculator.get_descriptor_primes()
print(temp_desc["H"][0].shape)
# print(temp_desc_prime["H"][0].shape)

descriptor_calculator.calculate_PCA(n_components = 10)

temp_desc = descriptor_calculator.get_descriptors()
# temp_desc_prime = descriptor_calculator.get_descriptor_primes()
print(temp_desc["H"][0].shape)
# print(temp_desc_prime["H"][0].shape)

descriptor_calculator.calculate_scaling()