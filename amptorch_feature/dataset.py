from torch.utils.data import Dataset
from .descriptor_base import AMPTorchDescriptorBase
from ase import Atoms

class AMPTorchDataset(Dataset):
    def __init__(
        self,
        images,
        descriptor,
        # label,
        forcetraining = True,
        store_fingerprints = True,
        parallel = False
        cores = 1,
    ):
        self.images = images
        self.forcetraining = forcetraining
        self.store_primes = store_primes
        
        assert isinstance(descriptor, AMPTorchDescriptorBase)
        # assert isinstance(images, Atoms)

        if not descriptor.fingerprints_ready:
            descriptor.prepare_fingerprints()
        
        if forcetraining and not descriptor.fingerprint_primes_ready:
            descriptor.prepare_fingerprints()
        

    def __len__(self):
        return len(self.atom_images)

    def __getitem__(self, index):
        pass