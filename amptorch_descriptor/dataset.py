from torch.utils.data import Dataset
from .descriptor_base import AMPTorchDescriptorBase
from ase import Atoms

class AMPTorchDataset(Dataset):
    def __init__(
        self,
        trajs,
        descriptor,
        # label,
        forcetraining = True,
        calculate_fingerprints = True,
        store_fingerprints = True,
        parallel = False,
        cores = 1,
    ):
        # self.images = images
        self.trajs = trajs
        self.descriptor = descriptor
        self.forcetraining = forcetraining
        self.calculate_fingerprints = calculate_fingerprints
        self.store_fingerprints = store_fingerprints

        self.fingerprints_ready = False
        self.fingerprint_primes_ready = False
        
        assert isinstance(descriptor, AMPTorchDescriptorBase)
        # assert isinstance(images, Atoms)

        if calculate_fingerprints:
            self.prepare_fingerprint_data()

    def prepare_fingerprint_data(self):
        self.descriptor.prepare_fingerprints(self.trajs, parallel=None, log=None, calculate_derivatives=self.forcetraining, save=self.store_fingerprints)
        self.fingerprints_ready = True
        if self.forcetraining:
            self.fingerprint_primes_ready = True

    def __len__(self):
        return len(self.atom_images)

    def __getitem__(self, index):
        pass