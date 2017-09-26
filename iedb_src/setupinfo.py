
import os

class SetupInfo(object):
    def __init__(self, version='20130222'):
        self.version = version
        self.path_main = '/fh/fast/gilbert_p/agartlan/gitrepo/iedb/mhc_i'
        self.path_method = os.path.join(self.path_main, 'method')
        self.path_data_base = os.path.join(self.path_main, 'data')
        self.path_data = os.path.join(self.path_data_base, 'MHCI_mhcibinding'+self.version) # This determines which version of training data to read.


