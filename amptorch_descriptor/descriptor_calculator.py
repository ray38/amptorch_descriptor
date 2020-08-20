# from torch.utils.data import Dataset
from .descriptor_base import AMPTorchDescriptorBase
from ase import Atoms
import pickle
from scipy.sparse import coo_matrix
import numpy as np

class DescriptorCalculator:
    def __init__(
        self,
        trajs,
        descriptor,
        automatic_calculation = True,
        calculate_descriptor_primes = True,
        store_descriptors = True,
        training_data = False,
        parallel = False,
        cores = 1,
    ):
        # self.images = images
        self.trajs = trajs
        self.descriptor = descriptor
        self.calculate_descriptor_primes = calculate_descriptor_primes
        self.training_data = training_data
        # self.calculate_fingerprints = calculate_fingerprints
        self.store_descriptors = store_descriptors

        self.element_list = self.descriptor._get_element_list()

        self.result_dir = "./_results_/"

        self.descriptors_ready = False
        # self.fingerprint_primes_ready = False
        
        assert isinstance(descriptor, AMPTorchDescriptorBase)
        # assert isinstance(images, Atoms)

        if automatic_calculation:
            self.prepare_descriptors()

    def prepare_descriptors(self):
        self.calculated_decsriptor_list = \
            self.descriptor.prepare_fingerprints(self.trajs, parallel=None, log=None, calculate_derivatives=self.calculate_descriptor_primes, save=self.store_descriptors)

        self.descriptors_ready = True

    def get_descriptors(self, separate_atomtypes = True):
        if self.descriptors_ready == False:
            print("ERROR, descriptors not calculated yet, please call prepare_descriptors() function first")
            return None
        
        if separate_atomtypes:
            result = {}
            for element in self.element_list:
                # result[element] = []
                element_descriptor_list = []
                for calculated_decsriptor in self.calculated_decsriptor_list:
                    temp = calculated_decsriptor[element]["descriptors"].copy()
                    element_descriptor_list.append(temp)
                    # element_descriptor_array = np.vstack([element_descriptor_array, temp]) if element_descriptor_array.size else temp
                result[element] = element_descriptor_array
            return result
        
        else: 
            print("WARNING: atomtype separation turned off, please make sure the dimensions match and you know what you are doing")
            result = []
            for calculated_decsriptor in self.calculated_decsriptor_list:
                descriptors = np.array([])
                for element in self.element_list:
                    temp = calculated_decsriptor[element]["descriptors"].copy()
                    descriptors = np.vstack([descriptors, temp]) if descriptors.size else temp
                result.append(descriptors)
            return result

    def get_descriptor_primes(self):
        pass


    def calculate_PCA(self, separate_atomtypes = True, save_models = True, n_components = 10, apply_PCA = True):
        from sklearn.decomposition import PCA

        if separate_atomtypes:
            models = {}
            raw_data = self.get_descriptors(separate_atomtypes=separate_atomtypes)
            for element in self.element_list:
                data = np.vstack(raw_data[element])
                pca_model = PCA(n_components=n_components)
                pca_model.fit(data)
                print(pca.explained_variance_ratio_)

                models[element] = pca_model
            if save:
                if not os.path.exists(self.desc_fp_database_dir):
                    os.makedirs(self.desc_fp_database_dir)
                pickle.dump( models, open( self.result_dir + "pca_models.p", "wb" ) )

            if apply_PCA:
                self.apply_PCA(models, separate_atomtypes=True)

        else:
            print("WARNING: NOT implemented")

    def apply_PCA(self, pca_model, separate_atomtypes = True):
        from sklearn.decomposition import PCA

        if separate_atomtypes:
            for element in self.element_list:
                model = pca_model[element]
                for calculated_decsriptor in calculated_decsriptor_list:
                    size_info = calculated_decsriptor[element]["size_info"]
                    calculated_decsriptor[element]["descriptors"] = model.transform(calculated_decsriptor[element]["descriptors"])
                    if calculate_descriptor_primes:
                        calculated_decsriptor[element]["descriptor_primes"] = \
                            self._apply_pca_model_to_descriptor_primes(calculated_decsriptor[element]["descriptor_primes"], size_info[0], size_info[1], size_info[2], model)
                    new_size_info = np.array([size_info[0], size_info[1], model.n_components_])
                    calculated_decsriptor[element]["size_info"] = new_size_info

        else:
            print("WARNING: NOT implemented")

    def _apply_pca_model_to_descriptor_primes(self, data, num_total_atoms, num_selected_atoms, num_original_descriptors, model):
        n_components = model.n_components_
        components = model.components_ #  shape (n_components, n_features)

        original_value = data["value"]
        original_row   = data["row"]
        original_col   = data["col"]
        original_size  = data["size"]

        # original_matrix = coo_matrix((data, (row, col)), shape=(4, 4)).toarray()
        original_array = coo_matrix((original_value, (original_row, original_col)), shape=original_size).toarray()
        print("original_size: {}".format(original_array.shape))

        transformed_array = np.zeros((num_selected_atoms * n_components, 3 * num_total_atoms))

        for i in range(num_selected_atoms):
            for j in range(n_components):
                index_new = i * n_components + j
                for k in range(num_original_descriptors):
                    index_old = i * num_original_descriptors + k
                    transformed_array[index_new] += components[j, k] * original_array[index_old]

        scipy_sparse_transformed = coo_matrix(transformed_array)
        print("density: {}%".format(100*len(scipy_sparse_fp_prime.data) / (fp_prime.shape[0] * fp_prime.shape[1])))
        transformed_result = {}
        transformed_result["value"] = scipy_sparse_transformed.data
        transformed_result["row"]   = scipy_sparse_transformed.row
        transformed_result["col"]   = scipy_sparse_transformed.col
        transformed_result["size"]  = np.array(transformed_array.shape)

        return transformed_result

    # def calculate_scale(self):

    # def apply_scale(self):

    # def get_descriptors(self):

    # def get_descriptor_primes(self):

    # def __len__(self):
    #     return len(self.atom_images)

    # def __getitem__(self, index):
    #     pass