from abc import ABC, abstractmethod
import numpy as np
import h5py
import torch
import os
import time
from .util import get_traj_hash
 
class AMPTorchDescriptorBase(ABC):
 
    def __init__(self):
        super().__init__()
        self.fp_database = "./_saved_fingerprints_/"

        # To Be specified/calculated
        self.descriptor_type = "Default"
        self.descriptor_setup_hash = "Default"


    
    @abstractmethod
    def calculate_fingerprints(self, image, params_set, calculate_derivatives=True):
        # image is a single snapshot
        pass

    @abstractmethod
    def get_descriptor_setup_hash(self):
        #set self.descriptor_setup_hash
        pass

    @abstractmethod
    def save_descriptor_setup(self, filename):
        pass

    @abstractmethod
    def prepare_descriptor_parameters(self):
        # prepare self.params_set
        pass
    
    def prepare_fingerprints(self, trajs, parallel=None, log=None, calculate_derivatives=True, save=True):
        # params_set = self.params_set

        Total_Num_Trajs = len(trajs)

        trajs_fingerprint_list = []
        trajs_fingerprint_prime_list = []

        if save:
            self.setup_fingerprint_database()

            for traj in trajs:
                traj_start_time = time.time()
                traj_hash = get_traj_hash(traj)
                traj_db_filename = "{}/AmpFP-{}-{}.h5".format(self.desc_fp_database_dir, self.descriptor_type, traj_hash)

                with h5py.File(traj_db_filename,'a') as db:

                    Total_Num_Snapshots = len(traj)
                    for i, snapshot in enumerate(traj):
                        start_time = time.time()
                        try:
                            current_snapshot_grp = db[str(i)]
                        except:
                            current_snapshot_grp = db.create_group(str(i))
                        
                        for element in self.elements:
                            if element in snapshot.get_chemical_symbols():
                                try:
                                    current_element_grp = current_snapshot_grp[element]
                                except:
                                    current_element_grp = current_snapshot_grp.create_group(element)

                                try:
                                    fps = np.array(current_element_grp["fps"])
                                    fp_primes_val = np.array(current_element_grp["fp_primes_val"])
                                    fp_primes_row = np.array(current_element_grp["fp_primes_row"])
                                    fp_primes_col = np.array(current_element_grp["fp_primes_col"])
                                    fp_primes_size = np.array(current_element_grp["fp_primes_size"])
                                except: 
                                    fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = \
                                        self.calculate_fingerprints(snapshot, element, calculate_derivatives=calculate_derivatives)

                                    current_element_grp.create_dataset("fps", data=fps)
                                    current_element_grp.create_dataset("fp_primes_val", data=fp_primes_val)
                                    current_element_grp.create_dataset("fp_primes_row", data=fp_primes_row)
                                    current_element_grp.create_dataset("fp_primes_col", data=fp_primes_col)
                                    current_element_grp.create_dataset("fp_primes_size", data=fp_primes_size)

                                indices = np.vstack((fp_primes_row, fp_primes_col))
                                torch_indices = torch.LongTensor(indices)
                                torch_values = torch.FloatTensor(fp_primes_val)
                                fp_prims_torch_sparse = torch.sparse.FloatTensor(torch_indices, torch_values, torch.Size(fp_primes_size))

                                trajs_fingerprint_list.append(torch.from_numpy(fps))
                                trajs_fingerprint_prime_list.append(trajs_fingerprint_list)
                            else:
                                print("element not in current snapshot: {}".format(element))
                        
                        took_time = time.time() - start_time
                        print("finished snapshot {}/{}, took time: {}".format(i+1, Total_Num_Snapshots, took_time))
                
                print("finished traj, took {}".format(time.time() - traj_start_time))

        else:
            for traj in list_of_trajs:
                for i, snapshot in enumerate(traj):
                    
                    for element in elements:
                        if element in snapshot.get_chemical_symbols():
                            fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = \
                                self.calculate_fingerprints(self, snapshot, element, calculate_derivatives=calculate_derivatives)

                            indices = np.vstack((fp_primes_row, fp_primes_col))
                            fp_prims_torch_sparse = torch.sparse.FloatTensor(indices, fp_primes_val, torch.Size(fp_primes_size))

                            trajs_fingerprint_list.append(torch.from_numpy(fps))
                            trajs_fingerprint_prime_list.append(trajs_fingerprint_list)
                        else:
                            print("element not in current snapshot: {}".format(element))


        return trajs_fingerprint_list, trajs_fingerprint_prime_list


    # def save_fingerprints(self, fingerprints, fingerprint_primes):
    #     # save and load as sparse matrices
        

    #     pass
    
    # def load_fingerprints(self):
    #     pass

    

    def setup_fingerprint_database(self):
        self.get_descriptor_setup_hash()
        self.desc_type_database_dir = "{}/{}".format(self.fp_database, self.descriptor_type)

        self.desc_fp_database_dir = "{}/{}".format(self.desc_type_database_dir, self.descriptor_setup_hash)

        if not os.path.exists(self.fp_database):
            os.makedirs(self.fp_database)
        
        if not os.path.exists(self.desc_type_database_dir):
            os.makedirs(self.desc_type_database_dir)

        if not os.path.exists(self.desc_fp_database_dir):
            os.makedirs(self.desc_fp_database_dir)

        descriptor_setup_filename = "__descriptor_setup__.txt"
        descriptor_setup_path = "{}/{}".format(self.desc_fp_database_dir, descriptor_setup_filename)
        self.save_descriptor_setup(descriptor_setup_path)
        
