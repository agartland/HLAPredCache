import os
import cPickle, re

from setupinfo import * #@UnusedWildImport

class InputError(Exception):
    """Exception raised for errors in the input."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value
    
class UnexpectedInputError(Exception):
    """Exception raised for errors in the input."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class PredictorError(Exception):
    """Exception raised for errors in the ARBMatrix and SMMMatrix classes."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class Proteins:
    """Contains a list of protein sequences and names and several conversion functions."""
    def __init__(self):
        self.sequences=[]
        self.names=[]

    def __len__(self):
        return len(self.sequences)

    def add_protein(self,sequence, length=None, source=None, name=None):
        sequence=sequence.strip().upper()
        if sequence != 'None' and length and length[0] != 'None':
            min_length = min([int(i) for i in length])
            if len(sequence) < min_length: 
                if name != None:
                    raise InputError("Input sequence %s in the %s is too short." %(name, source))
                else: raise InputError("Input sequence in the %s is too short." %source)

        for amino_acid in sequence:
            if not amino_acid in "ACDEFGHIKLMNPQRSTVWY":
                raise InputError("Sequence: '%s' contains an invalid character: '%c' at position %d." %(sequence, amino_acid, sequence.find(amino_acid)))
                
        self.sequences.append(sequence)
        
        if name==None:
            name = "sequence %d" %(len(self.sequences))
        self.names.append(name)

    def extract_fasta(self, fasta, length=None, source=None):
        input_sequences = fasta.split(">")
        if len(input_sequences)<2:
            raise InputError("Invalid Fasta format: No '>' found.")
        for i in input_sequences[1:]:
            if(len(i) > 0):
                end_of_name = i.find("\n")
                if end_of_name == -1:
                        raise InputError("Invalid Fasta format: No Protein sequence found between two names.")
                name = i[:end_of_name]
                seq = i[end_of_name:].split()
                self.add_protein("".join(seq), length, source, name)

class InputData(object):
    def __init__(self, version, method, mhc, hla_seq, length, proteins, species=None, protein_mhc=None, freq=None, negatives=None, duplicates=None, pred_tool=None):
        '''Fully define all inputs needed to run MHC class I binding predictions.'''
        self.set_version(version)
        self.set_method(method)
        self.set_mhc(mhc)
        self.set_pseudo(hla_seq)
        self.set_peptide_length(length)
        self.input_protein = proteins  # This contains the peptides that are tested for binding to a given mhc.
        self.input_protein_mhc = protein_mhc # This contains the user-provided mhc sequence.
        self.species=species
        self.freq = freq    # boolean data type (true/false)
        self.negatives = negatives
        self.duplicates = duplicates
        self.tool = pred_tool
        
    
    def get_freq(self): return self.freq
    
    def set_method(self, method): self.method = method

    def set_version(self, version): self.version = version

    def get_version(self): return self.version

    def set_mhc(self, mhc):
        '''Selected from a list.'''
        self.mhc = mhc
    
    def set_pseudo(self, hla_seq):
        '''User supplied MHC protein.'''
        self.hla_seq = hla_seq

    def set_peptide_length(self, length):
        self.length = length
        
    def set_input_protein(self, seq_list):
        '''Assumes a list of strings.'''
        self.input_protein = None
        temp_protein = Proteins()
        for seq in seq_list:
            temp_protein.add_protein(seq)
        self.input_protein = temp_protein

    def set_input_protein_mhc(self, usermhc_seq):
        self.mhc = 'UserDefined'
        self.input_protein_mhc = None
        temp_protein = Proteins()
        temp_protein.add_protein(usermhc_seq)
        self.input_protein_mhc = temp_protein


#class PredictorSelection(object):
#    '''Provides a coherent and easy to use mechanism to (1) validate selected (MHC,length) and (2) return selections.'''
#        # I want to make the following section of the code better. 20090629.
#        # It is not easy to read. It is hard to modify the list for each predictor.
#        # The list of (MHC,length) is not human readable.
#        # I don't like the structure of 'tools'. Why not a dictionary? instead of nested structure?
#        # [] Isolate this code segment for later updates.
#
#
#        # [] Validates the selection of (MHC,length). Makes sure it is present before making predictions.
#        # [] Defines what 'All alleles' means.
#        # [] Defines what 'All lengths' means.
#        # [] Returns a list of (MHC, length) that meet the selection criterion.
#        # [] Selection criterion: (method, MHC, length)
#        #        method = 'ann','arb',...
#        #        MHC = 'All alleles', 'HLA A*0201', ...
#        #        length = 'All lengths', 8, 9, ...
#    def __init__(self, pred_method, species, mhc, length):
#        self.pred_method = pred_method
#        self.species = species
#        self.mhc = mhc
#        self.length = length
#
#    def get_tool_selection(self):
#        tool_selection = []
#        tools = read_available_model_list()
#
#        # Check that the selected method is valid.
#        dic_method = get_method_dic()
#        if not dic_method.has_key(self.pred_method):
#            raise UnexpectedInputError("Selected prediction method '%s' does not exist." % (self.pred_method))
#
#        method_index = dic_method[self.pred_method]
#        'The user selected an mhc sequence from a prefined list.'
#        for (species, mhcs) in tools:
#            if species == self.species:
#                for (mhc, lengths) in mhcs:
#                    if self.mhc=="All alleles" or self.mhc == mhc:
#                        for (length, methods) in lengths:
#                            if self.length == "All lengths" or int(self.length) == int(length):
#                                if method_index in methods:
#                                    tool_selection.append((mhc, int(length)))
#        if len(tool_selection) == 0:
#            raise (UnexpectedInputError("Couldn't find tools matching species='%s', allele='%s', length='%s', method='%s'." % (self.species, self.mhc, self.length, self.pred_method)))
#        return tool_selection


class PredictorSelectionB(object):
    '''Provides a coherent and easy to use mechanism to (1) validate selected (MHC,length) and (2) return selections.'''
    def __init__(self, setupinfo):
        self.version = setupinfo.version
        ms = MethodSet()
        self.method_list = ms.get_method_list(version=self.version)
        self.method_list.sort()
        self.path_data = setupinfo.path_data
        self.dic_model_list = self.get_dic_model_list(self.method_list)

    def set_method_list(self, method_list):
        self.method_list = method_list
        self.dic_model_list = self.get_dic_model_list(method_list)

    def get_method_list(self):
        return self.method_list

    def get_dic_model_list(self, method_list):
        dic_model_list = {}
        for method in method_list:
            fname = os.path.join(self.path_data, method, 'model_list.txt')
            f=open(fname,'r')
            lines=f.readlines()
            model_list = [line.split('\t')[0].strip() for line in lines]
            f.close()
            dic_model_list[method] = model_list
        return dic_model_list

    def get_tool_selection(self, pred_method, mhc, length, species, tool=None):
        self.pred_method = pred_method
        self.species = species
        self.mhc = mhc
        self.length = length  # The peptide length of a query.
        
        mhc_all = 'All alleles'  # Used to select all available MHCs.
        length_all = 'All lengths' # Used to select all available lengths.

        tool_selection = []
        if (self.pred_method not in self.method_list):
            raise UnexpectedInputError("Selected prediction method '%s' does not exist." % (self.pred_method))
        
        
        processing_netmhcpan = [8,9,10,11,12,13,14]
        
        model_list = self.dic_model_list[self.pred_method] # A list of (MHC, length)
        for mhc_length in model_list:
            (mhc,length) = get_standard_mhc_name(mhc_length)
            species = get_species(mhc) #added
            if (self.species == species):
                if ((self.mhc == mhc) or (self.mhc == mhc_all)):
                    if (str(self.length) == str(length)) or (str(self.length) == length_all): # When comparing lengths, they should be of difference
                        if pred_method == 'IEDB_recommended' and tool == 'processing':
                            if int(length) in processing_netmhcpan:
                                tool_selection.append((mhc, int(length)))
                        else: tool_selection.append((mhc, int(length)))
                            
        if len(tool_selection) == 0:
            raise (UnexpectedInputError("Could not find tools matching species='%s', allele='%s', length='%s', method='%s'." % (self.species, self.mhc, self.length, self.pred_method)))
        return self.clean(tool_selection)

    def clean(self, lst):
        new_lst=[]
        for item in lst:
            if item !=0:
                if item not in new_lst:
                    new_lst.append(item)
        return new_lst

    def get_model_list(self, method):
        return self.dic_model_list[method]

    def get_model_list_combined(self):
        '''Returns a list of models such that at least one of the method in method_list is available.'''
        # Output: model_name, method_list_subset
        model_list = []  # In the format of a list of HLA-A-0201-9, etc.
        key_list = self.dic_model_list.keys()
        for method in key_list:
            model_list_temp = self.dic_model_list[method]
            model_list = model_list + model_list_temp
        model_list = list(set(model_list))
        return model_list

    def is_model_available(self, method, mhc_length):
        v=False
        model_list = self.dic_model_list[method]
        v = mhc_length in model_list
        return v

    def get_available_methods(self, mhc_length, method_list):
        method_list_available = []
        for method in method_list:
            v = self.is_model_available(method, mhc_length)
            if (v == True):
                method_list_available.append(method)
        return method_list_available

### Routines relavant to predictors.
# There should be only one place where these methods are listed.
# Q: What is the index of a method?
# Q: What is the official name of the method? e.g. NetMHC of 'ann'.
# Q: What methods are available for the version 20090901B?
class MethodSet(object):
    '''Holds available methods for each version of tools.'''
    def __init__(self):
        self.dic_method_info = self.get_dic_method_info()    # Stores (1) method_index, and (2) method_name_official for each method.
        self.dic_method_list = self.get_dic_method_list() # Stores a list of methods for a given version.

    def get_dic_method_info(self):
        '''List all methods. These namings and methods indexes should be constant across versions.'''
        dic = {}
        dic['ann']                = (1, 'NetMHC')            # method_name, method_index, official_name
        dic['smm']                = (2, 'SMM')
        dic['arb']                = (3, 'ARB')
        dic['smmpmbec']           = (4, 'SMMPMBEC')          # SMM with PMBEC covariance matrix.
        dic['comblib_sidney2008'] = (5, 'ComblibSidney2008') # Combinatorial matrices: for now using Sidney et al. 19 scoring matrices.
        dic['consensus']          = (6, 'Consensus')
        dic['netmhcpan']          = (7, 'NetMHCpan')         # NetMHCpan-2.0 for MHC-I
        #dic['comblib_udaka2000']  = (8, 'ComblibUdaka2000')  # Combinatorial matrices: using 4 scoring matrices by udaka et al.
        dic['pickpocket']          = (9, 'PickPocket')
        dic['netmhccons']          = (10, 'NetMHCcons')
        return dic

    def get_dic_method_list(self):
        '''The order of the methods is important.'''
        dic={}
        dic['20130222'] = ['ann', 'comblib_sidney2008', 'consensus', 'IEDB_recommended', 'netmhcpan', 'smm', 'smmpmbec', 'pickpocket', 'netmhccons']
        dic['20090901B'] = ['ann', 'arb', 'comblib_sidney2008', 'consensus', 'IEDB_recommended', 'netmhcpan', 'smm', 'smmpmbec']
        dic['20090901'] = ['ann', 'smm', 'arb', 'smmpmbec', 'comblib_sidney2008', 'consensus']
        dic['20071227'] = ['ann', 'smm', 'arb']
        dic['20060101'] = ['ann', 'smm', 'arb']
        return dic

    def get_method_index(self, method):
        if method != 'IEDB_recommended':  #method is checked for IEDB_recommended in the web tools
            (index, method_name_official) = self.dic_method_info[method]
            return index
        
#         (index, method_name_official) = self.dic_method_info[method]
#         return index

    def get_method_name_official(self, method):
        (index, method_name_official) = self.dic_method_info[method]
        return method_name_official

    def get_method_list(self, version='20130222'):
        return self.dic_method_list[version]

    def get_version_list(self):
        v = ['20130222', '20090901B', '20090901', '20071227', '20060101']
        return v

    def get_version_older(self, version='20130222'):
        '''Returns one version older. If starting with the oldest version, returns ths current version.'''
        version_list = self.get_version_list()
        version_index = version_list.index(version)
        version_older_index = version_index + 1 - len(version_list)  # 0 + 1 - 3=-2  # older version of 20090901B is 20090901 is 20071227
        version_older = version_list[version_older_index]
        return version_older

def get_path_model(path_data, mhc, length):
    '''Used by ARB, SMM to read appropriate files containing trained models.'''
    model_name = mhc.replace('*','-').replace(' ','-').replace(':','') + '-' + str(length)
    path_model = os.path.join(path_data, model_name) + '.cpickle'  # HLA-A-0201-9
    return path_model

### MHC naming related routines.
# Converts from 'HLA-A-0201-10' to '(HLA A*0201, 10)'
# How to deal with mouse:
# How to deal with patr
# How to deal with mamu
# How to deal with human
def get_standard_mhc_name(mhc_temp):
    temp = mhc_temp.strip().split('-')
    length = temp[-1]
    mhc = '-'.join(temp[0:-1])
#    mhc = mhc.replace('HLA-A','HLA-A*')
#    mhc = mhc.replace('HLA-B','HLA-B*')
#    mhc = mhc.replace('HLA-C','HLA-C*')
#    mhc = mhc.replace('HLA-E','HLA-E*')
#    mhc = mhc.replace('HLA-G','HLA-G*')
#    mhc = mhc.replace('H-2-','H-2-')
#    mhc = mhc.replace('Mamu-A-','Mamu-A*')
#    mhc = mhc.replace('Mamu-A','Mamu-A*') 
#    mhc = mhc.replace('Mamu-A1','Mamu-A1*')   # changed >> 'Mamu-A1*'
#    mhc = mhc.replace('Mamu-A2','Mamu-A2*')
#    mhc = mhc.replace('Mamu-A3','Mamu-A3*')
#    mhc = mhc.replace('Mamu-A4','Mamu-A4*')
#    mhc = mhc.replace('Mamu-A5','Mamu-A5*')
#    mhc = mhc.replace('Mamu-A6','Mamu-A6*')
#    mhc = mhc.replace('Mamu-A7','Mamu-A7*')
#    mhc = mhc.replace('Mamu-AG','Mamu-AG*')
#    mhc = mhc.replace('Mamu-B','Mamu-B*')
#    mhc = mhc.replace('Patr-A','Patr-A*')     # changed >> 'Patr-A*'
#    mhc = mhc.replace('Patr-B','Patr-B*')
#    mhc = mhc.replace('Patr-C','Patr-C*')
#    mhc = mhc.replace('SLA-1','SLA-1*')    # have 2 changed
#    mhc = mhc.replace('SLA-2','SLA-2*')    # have 2 changed
#    mhc = mhc.replace('SLA-3','SLA-3*')    # have 2 changed
#    mhc = mhc.replace('SLA-6','SLA-6*')
#    mhc = mhc.replace('Gogo-B','Gogo-B*')
#    mhc = mhc.replace('BoLA-N','BoLA-N*')
#    mhc = mhc.replace('BoLA-NC1','BoLA-NC1*')
#    mhc = mhc.replace('BoLA-NC2','BoLA-NC2*')
#    mhc = mhc.replace('BoLA-NC3','BoLA-NC3*')
#    mhc = mhc.replace('BoLA-NC4','BoLA-NC4*')
#    mhc = mhc.replace('*-','*')
#    mhc = mhc.replace('**','*')
    return (mhc, length)

def get_standard_mhc_name_b(mhc):
    mhc = mhc.replace('HLA-A-','HLA A*')
    mhc = mhc.replace('HLA-B-','HLA B*')
    mhc = mhc.replace('HLA-C-','HLA C*')
    mhc = mhc.replace('HLA-E-','HLA E*')
    mhc = mhc.replace('HLA-G-','HLA G*')
    mhc = mhc.replace('H-2-','H-2 ')
    mhc = mhc.replace('Mamu-A-','Mamu A*')
    mhc = mhc.replace('Mamu-B-','Mamu B*')
    mhc = mhc.replace('Patr-A-','Patr A*')
    mhc = mhc.replace('Patr-B-','Patr B*')
    mhc = mhc.replace('Gogo-B-','Gogo B*')
    return mhc

def get_species_list():
    '''List all species that are used here.'''
    species_list = ['chimpanzee', 'gorilla', 'mouse', 'macaque', 'pig', 'human', 'cow']
    return species_list

def get_species(mhc):
    species = None
    if re.search('HLA.*', mhc):
        species = 'human'
    elif re.search('H-2.*', mhc):
        species = 'mouse'
    elif re.search('Patr.*', mhc):
        species = 'chimpanzee'
    elif re.search('Mamu.*', mhc):
        species = 'macaque'
    elif re.search('Gogo.*', mhc):
        species = 'gorilla'
    elif re.search('SLA.*', mhc):
        species = 'pig'
    elif re.search('BoLA.*', mhc):
        species = 'cow'
    elif re.search('RT.*', mhc):
        species = 'rat'
    return species
### Math related routines. There should be a genral module for this type of routines. Something like R statistical language.
mean = lambda x: sum(x)/float(len(x))

def is_even(num):
    value = num % 2
    if (value == 0):
        return True
    else:
        return False

def median(v_temp):
    v = [x for x in v_temp]
    v.sort()
    value_median = None
    if (is_even(len(v)) == True):
        index_a = len(v)/2 - 1
        index_b = index_a + 1
        value_median = (v[index_a] + v[index_b])/2.0
    else:
        index = len(v)/2
        value_median = v[index]
    return value_median

def rankdata(v):
    # Rank starts from '1'
    # (value, its rank)
    # Those having same values get the average of their ranks.
    index_list = range(1,len(v)+1)
    v_a = [(x, index) for (x,index) in zip(v,index_list)]
    v_a.sort() # This will sort based on 'x' instead of 'index'
    v_rank = [(index, x, rank) for ((x,index), rank) in zip(v_a, index_list)]
    # 'rank' values should be now averaged if corresponding values are the same.
    v_rank.sort()  # Return the items to the original order.
    dic={}
    for (index,x,rank) in v_rank:
        key = x
        if (dic.has_key(key) == False):
            dic[key] = [rank]
        else:
            temp = dic[key]
            temp.append(rank)
    rank_list = [mean(dic.get(x)) for (index,x,rank) in v_rank]  # This averages ranks if their corresponding values are the same.
    return rank_list

def read_available_model_list():
    '''For each method, a list of available (mhc,length) is provided.'''
    setupinfo = SetupInfo()
    fname =setupinfo.fname_tool_list
    f=open(fname,'r')
    model_list = cPickle.load(f)
    f.close()
    return model_list

def read_available_model_list_b(version='20130222'):
    '''Builds tools list containing available predictors for all methods.'''
    setupinfo = SetupInfo(version=version)

    ms = MethodSet() # holds a list of methods for each version.
    method_list = ms.get_method_list(version=version)

    # The following block of code should be done in the very beginning of initialization of tools.
    dic_model_list={}
    for method in method_list:
        fname = os.path.join(setupinfo.path_data, method,'model_list.txt')
        f=open(fname,'r'); lines=f.readlines(); f.close()
        model_list = [line.split('\t')[0].replace('_','-').strip() for line in lines] # Assumes now 'HLA-A-0201-9 form.
        dic_model_list[method] = model_list

    # Group based on species, mhc, lengths. Each length is associated with a list of available methods.
    length_list = [8,9,10,11,12,13,14]
    species_list = get_species_list()
    model_list_all = []
    [model_list_all.extend(dic_model_list[method]) for method in method_list]
    model_list_all_unique = list(set(model_list_all))
    mhc_list = ['-'.join(name.split('-')[0:-1]) for name in model_list_all_unique]
    mhc_list_unique = list(set(mhc_list))

    tools = []
    for species in species_list:
        mhc_data_list = []
        for mhc in mhc_list_unique:
            mhc_old = mhc
            mhc = get_standard_mhc_name_b(mhc)
            species_test = get_species(mhc)
            if (species == species_test):
                length_data_list = []
                for length in length_list:
                    modelname = mhc_old+'-'+str(length)  # Should be 'HLA-A-0201-9'
                    method_index_list = []
                    for method in method_list:
                        model_list = dic_model_list[method]
                        method_index = ms.get_method_index(method)
                        if modelname in model_list:
                            method_index_list+=[method_index]
                    if (len(method_index_list) > 0):
                        length_data_list += [(int(length), method_index_list)]
                mhc_data_list.append((mhc, length_data_list))
        if (len(mhc_data_list)>0):
            mhc_data_list.sort()
            tools.append((species, mhc_data_list))
    return tools

def get_mhc_list(method_index):
    tools = read_available_model_list_b()
    tool_selection = []
    for (species, mhcs) in tools:
        for (mhc, lengths) in mhcs:
            for (length, methods) in lengths:
                if method_index in methods:
                    tool_selection.append((species, mhc, int(length)))
    return tool_selection

def get_peptides(sequence, peptide_length):
    '''Returns a list containing peptides of appropriate size.'''
    peptide_list = []
    num = len(sequence.strip()) - peptide_length + 1
    for i in range(num):
        peptide = sequence[i:i+peptide_length]
        peptide_list.append(peptide)
    return peptide_list

def list_mhc():
    #'''Returns a dictionary: key = method; value = vector of avaiable tools (ex. (HLA A*0201, 9))'''
    '''List all available mhc molecules.'''
    tools_list = []
    tools_list_a = get_mhc_list(1)
    tools_list_b = get_mhc_list(2)
    tools_list_c = get_mhc_list(3)
    tools_list.extend(tools_list_a)
    tools_list.extend(tools_list_b)
    tools_list.extend(tools_list_c)
    temp = set(tools_list)
    tools_list = list(temp)
    tools_list.sort()

    print "MHC", "\t", "Peptide Length"
    for row in tools_list:
        (species, mhc, peptidelength) = row
        print mhc, "\t", peptidelength
    return tools_list

def shared_list_mhc():
    '''Returns only those (MHC, peptide_length) combinations for which all methods are available. '''
    tools_list_a = get_mhc_list(1)
    tools_list_b = get_mhc_list(2)
    tools_list_c = get_mhc_list(3)
    aset = set(tools_list_a)
    bset = set(tools_list_b)
    cset = set(tools_list_c)
    shared_list_mhc = list(aset.intersection(bset).intersection(cset))
    shared_list_mhc.sort()
    print "MHC", "\t", "Peptide Length"
    for row in shared_list_mhc:
        (species, mhc, peptidelength) = row
        print mhc, "\t", peptidelength
    return shared_list_mhc

def test_predictor_selection():
    setupinfo = SetupInfo() # assumes version 20090901B.
    path_data = setupinfo.path_data
    pred_method = 'ann'
    species = 'human'
    mhc = 'HLA A*0201'
    length = '9'

    method_list = ['ann','smm','comblib_sidney2008', 'netmhcpan']  # This is the predictors used to build 'consensus'
    psb = PredictorSelectionB(path_data, pred_method, species, mhc, length)
    model_list = psb.get_model_list_combined()
    print 'debug model_list ', len(model_list)
    content = []
    for name in model_list:
        method_list_available = psb.get_available_methods(name, method_list)
        line = name+'\t'+str(method_list_available)
        content.append(line)
    content.sort()
    for line in content:
        print line

if __name__ == '__main__':
    test_predictor_selection()

