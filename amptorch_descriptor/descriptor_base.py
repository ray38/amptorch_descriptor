from abc import ABC, abstractmethod
from .util import get_traj_hash
import h5py
import os
 
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

        trajs_fingerprint_list = []
        trajs_fingerprint_prime_list = []

        if save:
            self.setup_fingerprint_database()

            for traj in list_of_trajs:
                traj_hash = get_traj_hash(traj)
                traj_db_filename = "{}/AmpFP-{}-{}.h5".format(self.desc_fp_database_dir, self.descriptor_type, traj_hash)

                with h5py.File(traj_db_filename,'a') as db:
                    for i, snapshot in enumerate(traj):
                        try:
                            current_snapshot_grp = db[i]
                        except:
                            current_snapshot_grp = db.create_group(i)
                        
                        for element in self.elements:
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
                                    self.calculate_fingerprints(self, snapshot, element, calculate_derivatives=calculate_derivatives)

                                current_element_grp.create_dataset("fps", data=fps)
                                current_element_grp.create_dataset("fp_primes_val", data=fp_primes_val)
                                current_element_grp.create_dataset("fp_primes_row", data=fp_primes_row)
                                current_element_grp.create_dataset("fp_primes_col", data=fp_primes_col)
                                current_element_grp.create_dataset("fp_primes_size", data=fp_primes_size)

                            indices = np.vstack((fp_primes_row, fp_primes_col))
                            fp_prims_torch_sparse = torch.sparse.FloatTensor(indices, fp_primes_val, torch.Size(fp_primes_size))

                            trajs_fingerprint_list.append(torch.from_numpy(fps))
                            trajs_fingerprint_prime_list.append(trajs_fingerprint_list)
        
        else:
            for traj in list_of_trajs:
                for i, snapshot in enumerate(traj):
                    
                    for element in elements:
                        fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = \
                            self.calculate_fingerprints(self, snapshot, element, calculate_derivatives=calculate_derivatives)

                        indices = np.vstack((fp_primes_row, fp_primes_col))
                        fp_prims_torch_sparse = torch.sparse.FloatTensor(indices, fp_primes_val, torch.Size(fp_primes_size))

                        trajs_fingerprint_list.append(torch.from_numpy(fps))
                        trajs_fingerprint_prime_list.append(trajs_fingerprint_list)


        return trajs_fingerprint_list, trajs_fingerprint_prime_list


    # def save_fingerprints(self, fingerprints, fingerprint_primes):
    #     # save and load as sparse matrices
        

    #     pass
    
    # def load_fingerprints(self):
    #     pass

    

    def setup_fingerprint_database(self):
        get_descriptor_setup_hash()
        self.desc_fp_database_dir = "{}/{}".format(self.fp_database, self.descriptor_type)

        if not os.path.exists(self.fp_database):
            os.makedirs(self.fp_database)
        
        if not os.path.exists(self.desc_fp_database_dir):
            os.makedirs(self.desc_fp_database_dir)

        descriptor_setup_filename = "__descriptor_setup__.txt"
        descriptor_setup_path = "{}/{}".format(desc_fp_database_dir, descriptor_setup_filename)
        self.save_descriptor_setup(descriptor_setup_path)
        
