import hashlib
import numpy as np
from ase.io.trajectory import Trajectory
from .constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT

# class Cosine_cutoff(object):
#     """Cosine functional form suggested by Behler.

#     Parameters
#     ---------
#     Rc : float
#         Radius above which neighbor interactions are ignored.
#     """

#     def __init__(self, Rc):

#         self.Rc = Rc

#     def __call__(self, Rij):
#         """
#         Parameters
#         ----------
#         Rij : float
#             Distance between pair atoms.

#         Returns
#         -------
#         float
#             The value of the cutoff function.
#         """
#         if Rij > self.Rc:
#             return 0.
#         else:
#             return 0.5 * (np.cos(np.pi * Rij / self.Rc) + 1.)

#     def prime(self, Rij):
#         """Derivative (dfc_dRij) of the Cosine cutoff function with respect to Rij.

#         Parameters
#         ----------
#         Rij : float
#             Distance between pair atoms.

#         Returns
#         -------
#         float
#             The value of derivative of the cutoff function.
#         """
#         if Rij > self.Rc:
#             return 0.
#         else:
#             return -0.5 * np.pi / self.Rc * np.sin(np.pi * Rij / self.Rc)

#     def todict(self):
#         return {'name': 'Cosine',
#                 'kwargs': {'Rc': self.Rc}}

#     def __repr__(self):
#         return ('<Cosine cutoff with Rc=%.3f from amp.descriptor.cutoffs>'
#                 % self.Rc)

def _gen_2Darray_for_ffi(arr, ffi, cdata="double"):
    # Function to generate 2D pointer for cffi  
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)
    return arr_p

def get_traj_hash(traj):
    
    # assert isinstance(traj, Trajectory)
    string = ""
    atoms = traj[0]
    string += str(temp_atoms.pbc)
    try:
        flattened_cell = atoms.cell.array.flatten()
    except AttributeError:  # older ASE
        flattened_cell = atoms.cell.flatten()
    for number in flattened_cell:
        string += "%.15f" % number
    for number in atoms.get_atomic_numbers():
        string += "%3d" % number
    for number in atoms.get_positions().flatten():
        string += "%.15f" % number
    
    md5 = hashlib.md5(string.encode("utf-8"))
    hash_result = md5.hexdigest()

    return hash_result

def list_symbols_to_indices(list_of_symbols):
    list_indices = []
    for symbol in list_of_symbols:
        list_indices.append(ATOM_SYMBOL_TO_INDEX_DICT[symbol])
    return np.array(list_indices, dtype=np.intc)

def list_indices_to_symbols(list_of_indices):
    list_symbols = []
    for index in list_of_indices:
        list_symbols.append(ATOM_INDEX_TO_SYMBOL_DICT[index])
    return list_symbols