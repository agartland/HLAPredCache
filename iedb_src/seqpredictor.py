import os
import math, tempfile, bisect

from util import * #@UnusedWildImport

class PredictorSet(object):
    '''A set of predictors. An appropriate set of predictors should be loaded for each version.'''
    def __init__(self, setupinfo):
        self.setupinfo = setupinfo
        self.path_method = setupinfo.path_method
        self.path_data = setupinfo.path_data
        self.version = setupinfo.version
        #self.dic_predictor  = self.get_dic_predictor()

    def get_dic_predictor(self):
        '''Read only those predictors which current version support.'''
        path_method = self.path_method
        path_data = self.path_data
        ms = MethodSet()
        method_list = ms.get_method_list(self.version)

        dic = {}
        for method in method_list:
            dic[method] = self.get_predictor(method, path_method, path_data)                      # Experimental
        return dic

    def get_predictor(self, method, method_used=None):
        path_method = self.path_method
        path_data = self.path_data
        predictor = None
        predictor = {
            'ann':               ANNPredictor(path_method, path_data),                               # Must be included.
            'smm':               SMMMatrix(path_method, path_data, method_name='smm', method_used=method_used),               # Must be included.
            'smmpmbec':          SMMMatrix(path_method, path_data, method_name='smmpmbec', method_used=method_used),          # Must be included.
            'arb':               ARBMatrix(path_method, path_data),                                  # Must be included.
            'netmhcpan':         NetMHCpanPredictor(path_method, path_data, method_used),                         # Experimental
            'comblib_sidney2008':CombinatorialLibrary(path_method, path_data, 'comblib_sidney2008'), # Must be included.
            'comblib_udaka2000': CombinatorialLibrary(path_method, path_data, 'comblib_udaka2000'),  # Must be included.
            'consensus':         ConsensusPredictor(self.setupinfo),                                 # Experimental
            'pickpocket':        PickPocketPredictor(path_method, path_data, method_used),
            'netmhccons':        NetMHCconsPredictor(path_method, path_data, method_used)
        }[method]
        return predictor

    def get_method_name_list(self):
        key_list = self.dic_predictor.keys()
        key_list.sort()
        return key_list

    def get_predictor_list(self, method_list):
        predictor_list = []
        for method in method_list:
            if (self.dic_predictor.has_key(method) == True):
                predictor = self.dic_predictor[method]
                predictor_list.append(predictor)
        return predictor_list


class ARBMatrix(object):
    """Can load and save ARB matrices in multiple formats and use them to score sequences. """
    def __init__(self, path_method, path_data):
        # path dependent loading of files takes place in 'initialize'.
        self.slope= None
        self.intercept= None
        self.length= None
        self.mat={}

        self.path_method = os.path.join(path_method, 'arb')
        self.path_data = os.path.join(path_data, 'arb')

    def initialize(self, mhc, length):
        self.mhc = mhc
        self.length = length
        self.path_model = get_path_model(self.path_data, mhc, length)
        self.pickle_load(self.path_model)

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'IC50 (nM)'

    def predict_sequence(self,sequence,pred,peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        return scores

    def predict_peptide_list(self, peptide_list):
        '''Given a list of peptides, return corresponding list of predicted binding scores.'''
        scores=[]
        for peptide in peptide_list:
            score = 0.0
            for pos in range(self.length):
                amino_acid=peptide[pos]
                try:
                    score+=self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'.""" % (amino_acid, peptide))
            score/=-self.length
            score-=self.intercept
            score/=self.slope
            score=math.pow(10,score)
            if score < 0.0001:    # Cap predictable values
                score = 0.0001
            elif score > 1e6:
                score = 1e6
            scores.append(score)
        return (tuple(scores))

    def loadARBTrainOutput(self, infile):
        infile=infile.read()
        lines=infile.split("\n")
        self.mat.clear()
        self.length=len(lines[0].split())-1
        for line in lines[0:20]:
            entries = line.split()
            numbers = []
            for e in entries[1:]:
                numbers.append(math.log10(float(e)))
            if len(numbers)!=self.length:
                raise PredictorError("Invalid number of columns in ARB matrix: " + str(len(numbers)), " expected: " + str(self.length) + ".")
            self.mat[line[0]]=tuple(numbers)
        p = infile.find("SLOPE")
        self.slope = float(infile[p + 5:infile.find("\n",p+1)])
        p = infile.find("INTERCEPT")
        self.intercept = float(infile[p + 9:infile.find("\n",p+1)])

    def pickle_dump(self, file_name):
        fout=open(file_name,"wb")
        cPickle.dump(self.length, fout)
        cPickle.dump(self.mat,fout)
        cPickle.dump(self.slope,fout)
        cPickle.dump(self.intercept,fout)
        fout.close()

    def pickle_load(self, file_name):
        fin = open(file_name,"rb")
        self.length = cPickle.load(fin)
        self.mat = cPickle.load(fin)
        self.slope = cPickle.load(fin)
        self.intercept = cPickle.load(fin)
        fin.close()


class SMMMatrix:
    """Can load and save SMM matrices in multiple formats and use them to score sequences """
    def __init__(self, path_method, path_data, method_name='smm', method_used=None):
        self.offset = None
        self.length = None
        self.mat={}
        self.method_used = method_used
        self.path_method = os.path.join(path_method, method_name)
        self.path_data = os.path.join(path_data, method_name)

    def initialize(self, mhc, length):
        self.mhc = mhc
        self.length = length
        self.path_model = get_path_model(self.path_data, mhc, length)
        self.pickle_load(self.path_model)

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'IC50 (nM)'

    def predict_sequence(self,sequence,pred,peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        
        #get percentile scores
        args = ('smm', self.mhc.replace("*",""), self.length)
        ps = PercentileScore(os.path.dirname(self.path_data), 'consensus', args)
        percentile = ps.get_percentile_score(scores)
        return zip(scores, percentile)

    def predict_peptide_list(self, peptide_list):
        scores=[]
        for peptide in peptide_list:
            score=self.offset
            for pos in range(self.length):
                amino_acid=peptide[pos]
                try:
                    score+=self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'.""" % (amino_acid, peptide))
            score=math.pow(10,score)
            scores.append(score)
        return (tuple(scores))

    def load_text_file(self, infile):
        lines=infile.readlines()
        self.mat.clear()
        self.length=int(lines[0].split()[1])
        for line in lines[1:21]:
            entries = line.split()
            numbers = []
            for e in entries[1:]:
                numbers.append(float(e))
            if len(numbers)!=self.length:
                raise PredictorError("Invalid number of columns in SMM matrix: " + str(len(numbers)), " expected: " + str(self.length) + ".")
            self.mat[line[0]]=tuple(numbers)
        self.offset = float(lines[21])

    def save_text_file(self, outfile):
        outfile.write("NumCols:\t" + str(self.length) +"\n")
        for letter in sorted(self.keys()):
            outfile.write(letter)
            for val in self.mat[letter]:
                outfile.write("\t" + str(val))
            outfile.write("\n")
        outfile.write(str(self.offset))

    def pickle_dump(self, file_name):
        fout = open(file_name,"wb")
        cPickle.dump(self.length, fout)
        cPickle.dump(self.mat,fout)
        cPickle.dump(self.offset,fout)
        fout.close()

    def pickle_load(self, file_name):
        fin = open(file_name,"rb")
        self.length = cPickle.load(fin)
        self.mat = cPickle.load(fin)
        self.offset = cPickle.load(fin)
        fin.close()
        

class ANNPredictor:
    def __init__(self, path_method, path_data):
        '''Predictor for Artificial Neural Network (aka. ann).'''
        self.path_executable = os.path.join(path_method, 'netMHC-3.4', 'netMHC')
        self.path_data = os.path.join(path_data, 'ann')

    def initialize(self, mhc, length):
        self.mhc = mhc.replace("*","")
        self.length = length

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'IC50 (nM)'

    def predict_sequence(self,sequence,pred,peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        scores = self.predict_peptide_list(sequence,peptideList)
        
        #get percentile scores
        args = ('ann', self.mhc.replace("*",""), self.length)
        ps = PercentileScore(os.path.dirname(self.path_data), 'consensus', args)
        percentile = ps.get_percentile_score(scores)
        return zip(scores, percentile)

    def predict_peptide_list(self, peptide_list,peptideList):
        '''This routine can be directly called so that you do not make a file for each prediction.'''
        infile = tempfile.NamedTemporaryFile(prefix=self.path_data+'/', suffix='input')
        if not peptideList:
            infile.write(">test\n%s\n" % peptide_list)
            infile.seek(0)
            cmd = self.path_executable + ' -a ' + self.mhc + ' -l ' + str(self.length) + ' ' + infile.name 
        else:
            for p in peptide_list:
                infile.write('%s\n' % p)
            infile.seek(0)
            cmd = self.path_executable + ' -a ' + self.mhc + ' -l ' + str(self.length) + ' -p ' + infile.name 

        f = os.popen(cmd)
        lines = f.readlines()
        pid = f.close()
        ic50s=[]
        for line in lines:
            row = line.split()
            if 'Entry_1' in row: # for the list with 'Entry_1' element
                ic50 = float(row[3])
                ic50s.append(ic50)
        
        if pid != None:  # According to python doc, return status of '0' is actually 'None' using os.popen.
            msg = "ANNPredictor did not execute. Path: %s, input :%s, Process ID: %d." % (self.path_executable, infile.name, pid)
            raise PredictorError(msg)
        
        infile.close()  # This gets automatically deleted.
        return tuple(ic50s)
    

class NetMHCpanPredictor:
    # Once the base dir path has been specified, the specification of the execuatable should be automatic.
    # Make sure folowing two directories have "all write" privileges!!
    # (1) /home/life/projects/iedb_tools_development/tools/netMHCpan-2.0/scratch  # This is defined in the netmhcpan script that came with the package.
    # (2) tool directory for netmhc. This is where files are being written by apache.
    def __init__(self, path_method, path_data, method_used=None):
        self.path_executable = os.path.join(path_method, 'netMHCpan-2.8', 'netMHCpan')
        self.path_scratch = os.path.join(path_data, 'netmhcpan') # The user should have a write access to this directory.
        self.path_data = os.path.join(path_data) 
        self.method = method_used
        
    def initialize(self, mhc, length, hla_seq=None):
        self.mhc = self.get_mhc_name(mhc)
        self.length = length
        self.hla_seq = hla_seq
    
    def read_score_distributions(self):
        fname = os.path.join(self.path_data, 'netmhcpan', 'distribution_netmhcpan_2_8_bin.cpickle')
        f=open(fname,'r')
        dic = cPickle.load(f)
        f.close()
        return dic
    
    def get_mhc_name(self, temp_mhc):
        if re.search('H-2.*', temp_mhc[0:3]):
            mhc = temp_mhc
        else:
            if re.search("HLA.*", temp_mhc):
                mhc = temp_mhc.replace("*","")
            elif re.search("SLA.*", temp_mhc) or re.search("Mamu.*", temp_mhc) or re.search("BoLA.*", temp_mhc):
                mhc = temp_mhc.replace("*",":")
            else:
                mhc = temp_mhc.replace("*","")
        return mhc

    def isReal(self, num):
        try:
            float(num)
            return True
        except ValueError:
            return False
    
    def parse_netmhcpan(self, content):
        '''A working version of a parser.'''
        scores = []
        for lines in content:
            if 'PEPLIST' in lines: 
                data_list = lines.split()
                if data_list[0].isdigit():
                #TODO: is this a duplicate?
                    if re.search("USER_DEF", data_list[1]):
                        peptide = data_list[2]
                        binding_affinity = float(data_list[4])
                        IC50_score = math.pow(50000,(1-binding_affinity))
                        scores.append(IC50_score)
                    else:
                        IC50_score = float(data_list[5])
                        scores.append(IC50_score)
        return scores

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'IC50 (nM)'

    def predict_sequence(self,sequence,pred,peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        ic50scores = []
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        if self.mhc != 'User-defined':
            if self.method == 'IEDB_recommended':
            	
            	# TODO: avoid search/replace, retrieve from the database instead 
                if re.search("SLA.*", self.mhc) or re.search("BoLA.*", self.mhc):
                    self.mhc = self.mhc.replace(":","_")
                    
                key = ('netmhcpan', self.mhc, self.length)
                self.dic_score_distributions = self.read_score_distributions()
                score_distribution = self.dic_score_distributions[key]
                ic50scores.append(scores)
                scores_netmhcpan = []
                scores_percentile = [self.get_percentile_score(score, score_distribution) for score in scores]
                for score in scores_percentile:
                    scores_netmhcpan.append(score)
                    #raise PredictorError("this means you skipped the get_percentile_score")
                return tuple(scores_netmhcpan), ic50scores
            else: return scores
        else:
            return scores
    
    def search(self, a, x):
        'Find leftmost value greater than x'
        i = bisect.bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return a[i-1]
        #raise ValueError

    def get_percentile_score(self, score, score_distributions):
        '''For each score in score_list, what percentage of the scores in score_distributions is worse?
        Smaller the score, the more significant.'''
        score_percentile = None
        percentile = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
                      2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 
                      4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 
                      6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
                      8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 
                      11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
                      31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
                      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                      71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
                      91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
        right_dist_score = self.search(score_distributions, score)
        
        if score not in score_distributions:
            right_indx = score_distributions.index(right_dist_score)
            return percentile[right_indx]
        else:
            return percentile[score_distributions.index(score)]
    
    def predict_peptide_list(self, peptide_list):
        '''This routine can be directly called so that you do not make a file for each prediction.'''
        
        infile=tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='input')
        infile.close()
        infile = open(infile.name, "w")
        for peptide in peptide_list:
            infile.write(peptide + "\n")
        infile.close()
        
        usermhcfile = tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='usermhc')
        usermhcfile.close()
        
        if self.hla_seq == '': self.hla_seq = None
        if self.hla_seq != None:
            user_input = self.hla_seq.split(">")
            for i in user_input[1:]:
                if(len(i) > 0):
                    end_of_name = i.find("\n")
                    name = i[:end_of_name].rstrip()
                    name = '>' + name
                    seq = i[end_of_name:].split()
                    usermhcseq = "".join(seq)
                    usermhc = name + '\n' + usermhcseq
            usermhcfile = open(usermhcfile.name, "wb")
            usermhcfile.write(usermhc)
            usermhcfile.close()
        # Added by AFG on Feb 16, 2017: not sure whythis was needed today, but it works now
        tmpdirCMD = ' -tdir ' + os.path.join(self.path_scratch, 'tmp')        
        if self.hla_seq is None:
            cmd = self.path_executable + ' -p ' + infile.name + ' -ic50 -a ' + self.mhc + ' -l ' + str(self.length) + tmpdirCMD
        else:
            cmd = self.path_executable + ' -p ' + infile.name  + ' -hlaseq ' + usermhcfile.name + ' -l ' + str(self.length) + tmpdirCMD
        
        f = os.popen(cmd)
        
        content = f.readlines()
        scores = self.parse_netmhcpan(content)
        pid = f.close()
        
        os.remove(infile.name)
        if(self.hla_seq != None):
            os.remove(usermhcfile.name)
        
        if (len(peptide_list) != len(scores)):
            msg = ""
            for row in content:
                if (row[0] != "#"):
                    msg = msg + row + '<br>'
            raise PredictorError(msg)

        if cmd == None:  # According to python doc, return status of '0' is actually 'None' using os.popen.
            msg = "NetMHCpanPredictor did not execute. Path: %s, input :%s." % (self.path_executable, infile.name)
            raise PredictorError(msg)
        return (tuple(scores))

class PickPocketPredictor:
    
    def __init__(self, path_method, path_data, method_used=None):
        self.path_executable = os.path.join(path_method, 'pickpocket-1.1', 'PickPocket')
        self.path_scratch = os.path.join(path_data, 'pickpocket') # The user should have a write access to this directory.
        self.path_data = os.path.join(path_data) 
        self.method = method_used
        
    def initialize(self, mhc, length, hla_seq=None):
        self.mhc = self.get_mhc_name(mhc)
        self.length = length
        self.hla_seq = hla_seq
#     
    def get_mhc_name(self, temp_mhc):
        if re.search('H-2.*', temp_mhc[0:3]):
            mhc = temp_mhc
        else:
            if re.search("SLA.*", temp_mhc) or re.search("Mamu.*", temp_mhc):
                mhc = temp_mhc.replace("*",":")
            else:
                mhc = temp_mhc.replace("*","")
        return mhc
# 
    def transform_score(self, x):
        '''Given nethmcpan score, returns log10(IC50).
        1-log50k(aff) = score
        aff = 50000^(-(score - 1))'''
        ba = float(x)
        try:
            ic50 = math.pow(50000,(1-ba))
        except:
            print 'EXCEPTION, math error, x', x
        return ic50
#     

    def parse_pickpocket(self, content):
        '''A working version of a parser.'''
        scores = []
        for lines in content:
            if 'PEPLIST' in lines: 
                data_list = lines.split()
                if data_list[0].isdigit():
                #TODO: is this a duplicate?
                    if re.search("USER_DEF", data_list[1]):
                        peptide = data_list[2]
                        binding_affinity = float(data_list[4])
                        IC50_score = math.pow(50000,(1-binding_affinity))
                        scores.append(IC50_score)
                    else:
                        IC50_score = self.transform_score(data_list[4])
                        scores.append(IC50_score)
        return scores
# 

    def predict_sequence(self,sequence,pred, peptideList = False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        ic50scores = []
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        
        if self.mhc != 'User defined':
            if re.search("BoLA.*", self.mhc) or re.search("Mamu.*", self.mhc):
                args = ('pickpocket', self.mhc.replace(":","_"), self.length)
            elif re.search("SLA.*", self.mhc):
                key_mhc = re.sub(r'(SLA-\d)(:)(.+)', r'\1_\3', self.mhc)
                key_mhc = re.sub(r'(SLA-\d)(-)(.+)', r'\1\3', key_mhc)
                args = ('pickpocket', key_mhc, self.length)
            else:
                args = ('pickpocket', self.mhc.replace("*",""), self.length)
            ps = PercentileScore(self.path_data, 'pickpocket', args)
            percentile = ps.get_percentile_score(scores)
            return zip(scores, percentile)
        else:
            return scores
    
    def predict_peptide_list(self, peptide_list):
        '''This routine can be directly called so that you do not make a file for each prediction.'''
        infile=tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='input')
        infile.close()
        infile = open(infile.name, "w")
        for peptide in peptide_list:
            infile.write(peptide + "\n")
        infile.close()
 
        usermhcfile = tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='usermhc')
        usermhcfile.close()
         
        if self.hla_seq == '': self.hla_seq = None
        if self.hla_seq != None:
            user_input = self.hla_seq.split(">")
            for i in user_input[1:]:
                if(len(i) > 0):
                    end_of_name = i.find("\n")
                    name = i[:end_of_name].rstrip()
                    name = '>' + name
                    seq = i[end_of_name:].split()
                    usermhcseq = "".join(seq)
                    usermhc = name + '\n' + usermhcseq
            usermhcfile = open(usermhcfile.name, "wb")
            usermhcfile.write(usermhc)
            usermhcfile.close()
        
        if self.hla_seq is None:
            cmd = self.path_executable + ' -a ' + self.mhc + ' -l ' + str(self.length) + ' -inptype 1 -p ' + infile.name 
        else:
            cmd = self.path_executable + ' -hlaseq ' + usermhcfile.name + ' -l ' + str(self.length) + ' -inptype 1 -p ' + infile.name 
#             cmd = self.path_executable + ' -p ' + infile.name  + ' -hlaseq ' + usermhcfile.name + ' -l ' + str(self.length)
         
        f = os.popen(cmd)
        content = f.readlines()
        scores = self.parse_pickpocket(content)
        pid = f.close()
         
        os.remove(infile.name)
        if(self.hla_seq != None):
            os.remove(usermhcfile.name)
         
        if (len(peptide_list) != len(scores)):
            msg = ""
            for row in content:
                if (row[0] != "#"):
                    msg = msg + row + '<br>'
            raise PredictorError(msg)
 
        if cmd == None:  # According to python doc, return status of '0' is actually 'None' using os.popen.
            msg = "NetMHCpanPredictor did not execute. Path: %s, input :%s." % (self.path_executable, infile.name)
            raise PredictorError(msg)
        return (tuple(scores))
    

class NetMHCconsPredictor:
    
    def __init__(self, path_method, path_data, method_used=None):
        self.path_executable = os.path.join(path_method, 'netMHCcons-1.1', 'netMHCcons')
        self.path_scratch = os.path.join(path_data, 'netmhccons') # The user should have a write access to this directory.
        self.path_data = os.path.join(path_data) 
        self.method = method_used
        
    def initialize(self, mhc, length, hla_seq=None):
        self.mhc = self.get_mhc_name(mhc)
        self.length = length
        self.hla_seq = hla_seq

    def get_mhc_name(self, temp_mhc):
        if re.search('H-2.*', temp_mhc[0:3]):
            mhc = temp_mhc
        else:
            if re.search("SLA.*", temp_mhc) or re.search("Mamu.*", temp_mhc):
                mhc = temp_mhc.replace("*",":")
            else:
                mhc = temp_mhc.replace("*","")
        return mhc

    def parse_pickpocket(self, content):
        '''A working version of a parser.'''
        scores = []
        for lines in content:
            if 'PEPLIST' in lines: 
                data_list = lines.split()
                if data_list[0].isdigit():
                    scores.append(float(data_list[5]))
        return scores

    def predict_sequence(self,sequence,pred, peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        ic50scores = []
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        
        if self.mhc != 'User defined':
            if re.search("BoLA.*", self.mhc) or re.search("Mamu.*", self.mhc):
                args = ('pickpocket', self.mhc.replace(":","_"), self.length)
            elif re.search("SLA.*", self.mhc):
                key_mhc = re.sub(r'(SLA-\d)(:)(.+)', r'\1_\3', self.mhc)
                key_mhc = re.sub(r'(SLA-\d)(-)(.+)', r'\1\3', key_mhc)
                args = ('pickpocket', key_mhc, self.length)
            else:
                args = ('pickpocket', self.mhc.replace("*",""), self.length)
            ps = PercentileScore(self.path_data, 'pickpocket', args)
            percentile = ps.get_percentile_score(scores)
            return zip(scores, percentile)
        else:
            return scores
        
    
    def predict_peptide_list(self, peptide_list):
        '''This routine can be directly called so that you do not make a file for each prediction.'''
        infile=tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='input')
        infile.close()
        infile = open(infile.name, "w")
        for peptide in peptide_list:
            infile.write(peptide + "\n")
        infile.close()
 
        usermhcfile = tempfile.NamedTemporaryFile(prefix=self.path_scratch, suffix='usermhc')
        usermhcfile.close()
         
        if self.hla_seq == '': self.hla_seq = None
        if self.hla_seq != None:
            user_input = self.hla_seq.split(">")
            for i in user_input[1:]:
                if(len(i) > 0):
                    end_of_name = i.find("\n")
                    name = i[:end_of_name].rstrip()
                    name = '>' + name
                    seq = i[end_of_name:].split()
                    usermhcseq = "".join(seq)
                    usermhc = name + '\n' + usermhcseq
            usermhcfile = open(usermhcfile.name, "wb")
            usermhcfile.write(usermhc)
            usermhcfile.close()
             
        if self.hla_seq is None:
            cmd = self.path_executable + ' -a ' + self.mhc + ' -length ' + str(self.length) + ' -inptype 1 -f ' + infile.name 
        else:
            cmd = self.path_executable + ' -hlaseq ' + usermhcfile.name + ' -length ' + str(self.length) + ' -inptype 1 -f ' + infile.name 
         
        f = os.popen(cmd)
        content = f.readlines()
        scores = self.parse_pickpocket(content)
        pid = f.close()
         
        os.remove(infile.name)
        if(self.hla_seq != None):
            os.remove(usermhcfile.name)
         
        if (len(peptide_list) != len(scores)):
            msg = ""
            for row in content:
                if (row[0] != "#"):
                    msg = msg + row + '<br>'
            raise PredictorError(msg)
 
        if cmd == None:  # According to python doc, return status of '0' is actually 'None' using os.popen.
            msg = "netMHCconsPredictor did not execute. Path: %s, input :%s." % (self.path_executable, infile.name)
            raise PredictorError(msg)
        return (tuple(scores))
    
class CombinatorialLibrary:
    """Can load and save SMM matrices in multiple formats and use them to score sequences """
    def __init__(self, path_method, path_data, lib_source):
        self.dic_pssm = None
        self.offset = None
        self.length = None
        self.mat={}

        self.path_method = path_method
        self.path_data = path_data
        self.lib_source = lib_source

    def initialize(self, mhc, length):
        self.dic_pssm = self.read_pssm_comblib(self.lib_source)
        self.mhc = mhc
        self.length = length
        if re.search('H-2.*', self.mhc):
            i = re.search('(?<=\d)', self.mhc).start()
            key = (self.mhc[:i]+self.mhc[i:].replace("-","_"), self.length)
        else:
            key = (self.mhc.replace('-','_').replace('*','-').replace(':',''), self.length)
        
        w = self.dic_pssm[key]
        (self.mat, self.offset) = self.get_dic_mat(w)

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'Score'

    def predict_sequence(self,sequence,pred,peptideList=False):
        '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
        if not peptideList:
            peptide_list = get_peptides(sequence, self.length)
        else:
            peptide_list = sequence
        scores = self.predict_peptide_list(peptide_list)
        
        #get percentile scores
        args = ('comblib_sidney2008', self.mhc.replace("*",""), self.length)
        ps = PercentileScore(self.path_data, 'consensus', args)
        percentile = ps.get_percentile_score(scores)
        return zip(scores, percentile)

    def predict_peptide_list(self, peptide_list):
        scores = []
        for peptide in peptide_list:
            score = self.offset
            for pos in range(self.length):
                amino_acid = peptide[pos]
                try:
                    score += self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'.""" % (amino_acid, peptide))
            score = math.pow(10,score)
            scores.append(score)
        return (tuple(scores))

    def read_data_cpickle(self,fname):
        f = open(fname,'r')
        data = cPickle.load(f)
        f.close()
        return data

    def read_pssm_comblib(self, lib_source):
        'Reads in all available pssms derived from combinatorial libraries.'
        factor = 1.0 # This will be multipled to all matrix elements.
        fname_sidney2008 = os.path.join(self.path_data,'comblib_sidney2008','dic_pssm_sidney2008.cPickle')
#         fname_udaka2000 = os.path.join(self.path_data,'comblib_udaka2000','dic_pssm_udaka2000.cPickle')
        dic_pssm_sidney2008 = self.read_data_cpickle(fname_sidney2008)
#         dic_pssm_udaka2000 = self.read_data_cpickle(fname_udaka2000)
        dic_pssm = None
        if (lib_source == 'comblib_sidney2008'):
            factor = -1.0
            dic_pssm = dic_pssm_sidney2008
#         elif (lib_source == 'comblib_udaka2000'):
#             factor = 1.0
#             dic_pssm = dic_pssm_udaka2000

        key_list = dic_pssm.keys()
        for key in key_list:
            w = dic_pssm[key]
            w = [factor*val for val in w]
            dic_pssm[key] = w
        return dic_pssm

    def get_dic_mat(self, w):
        'Converts 1-dimensional vector into a dictionary of lists key = [aa]'
        offset = w[0]
        dic_mat = {}
        aa_list = "ACDEFGHIKLMNPQRSTVWY"
        for aa_index in range(len(aa_list)):
            aa = aa_list[aa_index]
            row = []
            for pos_index in range(self.length):
                index = 1 + 20*pos_index + aa_index
                value = w[index]
                row.append(value)
            dic_mat[aa] = row
        return (dic_mat, offset)

    def pickle_dump(self, file_name):
        fout = open(file_name,"wb")
        cPickle.dump(self.length, fout)
        cPickle.dump(self.mat,fout)
        cPickle.dump(self.offset,fout)
        fout.close()

    def pickle_load(self, file_name):
        fin = open(file_name,"rb")
        self.length = cPickle.load(fin)
        self.mat = cPickle.load(fin)
        self.offset = cPickle.load(fin)
        fin.close()


class ConsensusPredictor(object):
    '''Should consensus return only its scores, or those of other predictors as well?'''
    def __init__(self, setupinfo):
        self.setupinfo = setupinfo
        self.path_method = setupinfo.path_method
        self.path_data = setupinfo.path_data
        self.method_name_list = ['ann', 'smm', 'comblib_sidney2008'] #'netmhcpan' # Include these tools for each (MHC,length) if possible.
        self.score_array = [None, None] # (scores_predictor, scores_consensus) Holds predictor specific prediction scores.

    def initialize(self, mhc, length):
        self.dic_score_distributions = self.read_score_distributions()
        self.dic_predictor = self.get_dic_predictor()
        self.predictor_selection = PredictorSelectionB(self.setupinfo)
        self.predictor_selection.set_method_list(self.method_name_list) # This limits the list of methods available to those that consecern consensus.
        
        self.mhc = mhc
        self.length = length
        mhc_length = mhc + '-' + str(length)
        
        # Q: What methods are available for (mhc,length)?
        self.available_method_list = self.predictor_selection.get_available_methods(mhc_length, self.method_name_list)
        self.predictor_list = [self.dic_predictor[method] for method in self.available_method_list] # Get only those predictors for which (mhc,length) is available.
        # [] Build a list of predictors based on this set.
        [predictor.initialize(mhc,length) for predictor in self.predictor_list]

    def read_score_distributions(self):
        fname = os.path.join(self.path_data, 'consensus', 'distribution_consensus_bin.cpickle')
        f=open(fname,'r')
        dic = cPickle.load(f)
        f.close()
        return dic

    def get_score_array(self):
        return self.score_array

    def get_score_unit(self):
        '''The unit of prediction scores'''
        return 'Percentile'

    def predict_sequence(self,sequence,pred,peptideList=False):
        scores_predictor = [] # Lower the score, the better.
        ic50scores = []
        for (predictor, method_name) in zip(self.predictor_list, self.available_method_list):
            self.mhc = self.mhc.replace("*","")
            key = (method_name, self.mhc, self.length)
            score_distribution = self.dic_score_distributions[key]
            
            if method_name == 'smm' or method_name == 'ann' or method_name == 'comblib_sidney2008':
                spercentile = predictor.predict_sequence(sequence,pred,peptideList)
                scores = tuple([sp[0] for sp in spercentile])
                percentile = tuple([sp[1] for sp in spercentile])
                ic50scores.append(scores)
                scores_predictor.append(percentile)
            else:
                scores = predictor.predict_sequence(sequence,pred,peptideList)
                ic50scores.append(scores)
                # here scores = individual scores for each of the methods
                scores_percentile = [self.get_percentile_score(score, score_distribution) for score in scores]  #range = [0....100]
                scores_predictor.append(tuple(scores_percentile))
        
        scores_consensus = []
        for i in range(len(scores_predictor[0])):
            score_row = [scores[i] for scores in scores_predictor]
            scores_consensus.append(median(score_row))
            
        self.score_array = (scores_predictor, scores_consensus)
        if pred == 'submit_processing':
            return tuple(scores_consensus)
        else:
            ic50_ranks = zip(ic50scores, scores_predictor)
            return tuple(scores_consensus), ic50_ranks    #ic50scores
        
    def search(self, a, x):
        'Find leftmost value greater than x'
        i = bisect.bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return a[i-1]
        #raise ValueError
          
    def get_percentile_score(self, score, score_distributions):
        '''For each score in score_list, what percentage of the scores in score_distributions is worse?
        Smaller the score, the more significant.'''
        score_percentile = None
        percentile = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
                      2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 
                      4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 
                      6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
                      8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 
                      11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
                      31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
                      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                      71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
                      91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
        right_dist_score = self.search(score_distributions, score)
        
        if score not in score_distributions:
            right_indx = score_distributions.index(right_dist_score)
            return percentile[right_indx]
        else:
            return percentile[score_distributions.index(score)]

    def get_dic_predictor(self):
        path_method = self.path_method
        path_data = self.path_data
        dic = {}
        dic['ann']                = ANNPredictor(path_method, path_data) # Must be included.
        dic['smm']                = SMMMatrix(path_method, path_data) # Must be included.
        dic['arb']                = ARBMatrix(path_method, path_data) # Must be included.
        dic['netmhcpan']          = NetMHCpanPredictor(path_method, path_data) # Experimental
        dic['comblib_sidney2008'] = CombinatorialLibrary(path_method, path_data, 'comblib_sidney2008')
        dic['comblib_udaka2000']  = CombinatorialLibrary(path_method, path_data, 'comblib_udaka2000')
        #dic['consensus']          = ConsensusPredictor(path_method, path_data) # Experimental
        return dic


class PercentileScore:
    '''Should consensus return only its scores, or those of other predictors as well?'''
    def __init__(self, path_data, method, args):
        self.path_data = path_data
        self.method = method
        self.key = args
        self.score_distributions = self.read_score_distributions()[self.key]
        
    
    def get_percentile_score(self, scores):
        return tuple([self.scores(score) for score in scores])  #range = [0....100]
    
    
    def scores(self, score):
        '''For each score in score_list, what percentage of the scores in score_distributions is worse?
        Smaller the score, the more significant.'''
        score_percentile = None
        # TODO: this list and the one for consensus method should point to a same single list (or query from db?)
        percentile = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
                      2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 
                      4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 
                      6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
                      8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 
                      11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
                      31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
                      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                      71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
                      91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
        right_dist_score = self.search(self.score_distributions, score)
        
        if score not in self.score_distributions:
            right_indx = self.score_distributions.index(right_dist_score)
            return percentile[right_indx]
        else:
            return percentile[self.score_distributions.index(score)]
    
    def search(self, a, x):
        'Find leftmost value greater than x'
        i = bisect.bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return a[i-1]
        #raise ValueError
        
    def read_score_distributions(self):
        _fname = 'netmhcpan_2_8' if self.method == 'netmhcpan' else self.method
        fname = os.path.join(self.path_data, self.method, 'distribution_{0}_bin.cpickle'.format(_fname))
        f=open(fname,'r')
        dic = cPickle.load(f)
        f.close()
        return dic


class MHCBindingPredictions:
    '''A main class to access different MHC binding prediction methods.
       Also provides various types of checks such as whether a selected mhc molecule is accessible.'''
    def __init__(self, input):
        #== The following set of variables should be cleaned up.
        self.version       = input.version
        self.method        = input.method    # User-selected predictive method.
        self.method_set_selected = []        # To be used when 'iedb recommended' used; methods used across a list of (mhc,length).
        self.species       = input.species
        self.length        = input.length
        self.proteins      = input.input_protein        # source of peptides to make binding predictions.
        self.proteins_mhc  = input.input_protein_mhc    # User-provided mhc sequence.
        self.setupinfo     = SetupInfo(version=self.version)
        self.path_data     = self.setupinfo.path_data
        self.predictor_set = PredictorSet(self.setupinfo)
        self.mhc           = input.mhc   # mhc allele
        self.hla_seq       = input.hla_seq   # user input mhc sequence
        self.freq          = input.freq  # boolean data type (true/false)
        self.negatives      = input.negatives  # boolean data type (true/false)
        self.duplicates    = input.duplicates  # list of duplicate allele-length pairs
        self.tool          = input.tool
        
        if(self.mhc == 'Allele'):
            self.mhc = None
        
        ps = PredictorSelectionB(self.setupinfo)
        self.tool_selection = []
        if self.hla_seq == '': self.hla_seq = None
        if self.hla_seq is None:
            for m, l, s in zip(self.mhc, self.length, self.species):
                self.tool_selection.extend(ps.get_tool_selection(self.method, m, l, s, self.tool))

    def get_score_unit(self):
        '''The unit of prediction scores'''
        prediction_score_unit = 'ic50'
        if (self.method=='comblib_sidney2008'):
            prediction_score_unit = 'score'
        elif self.method=='consensus':
            prediction_score_unit = 'consensus_percentile_rank'
        elif self.method=='IEDB_recommended':
            prediction_score_unit = 'percentile_rank'
        return prediction_score_unit
    
    def get_method_set_selected(self, method):
        method_set_selected = []
        if method == 'IEDB_recommended' or method == 'consensus':
            for allele, length in self.tool_selection:
                method_set_selected.extend(self.method_lookup(allele, length))
            method_set_selected = list(set(method_set_selected))
            num_methods = len(method_set_selected)
            if num_methods >= 2: method_set_selected.append('consensus')
            method_set_selected.sort()
        else: method_set_selected.append(method)
        return method_set_selected

    def predict(self, sequence_list, pred=None, peptideList = False):
        results = []
        
        if self.method == 'IEDB_recommended':
            for allele, length in self.tool_selection:
                method_name = ','.join(self.method_lookup(allele, length))
                # if 'recommended' is passed from MHC-I processing, overwrite the method_name to 'netmhcpan'
                if pred == 'submit_processing':
                    method_name = 'netmhcpan'
                     
                predictor = ''
                if method_name != 'netmhcpan':
                    predictor = self.predictor_set.get_predictor('consensus', 'IEDB_recommended')
                else:
                    predictor = self.predictor_set.get_predictor('netmhcpan', 'IEDB_recommended')
                     
                predictor.initialize(allele, int(length))
                scores = []
                if not peptideList:
                    for sequence in sequence_list:
                        scores.append(predictor.predict_sequence(sequence,pred))
                else:
                    scores.append(predictor.predict_sequence(sequence_list,pred,peptideList))
                results.append((length, allele, scores, method_name))
        else:
            
            #if 'IEDB recommended' is not selected, only one method will be chosen:
            predictor = self.predictor_set.get_predictor(self.method)
            
            self.method_set_selected = [self.method]
            if ((self.hla_seq is None) | (self.hla_seq == "")):
                '''
                1. A predictive method was chosen.
                2. Loop over a set of (allele,length) alleles selected by the user.
                3. For each combination: 
                    Loop over each sequence and make predictions.
                4. Collect scores as a list.
                '''
                for allele, length in self.tool_selection:
                    method_name = ','.join(self.method_lookup(allele, length))
                 
                    if not method_name:
                        method_name = 'netmhcpan'
 
                    if (self.method == 'netmhcpan'): predictor.initialize(allele, int(length), self.hla_seq)
                    else: predictor.initialize(allele, int(length))
 
                    scores = []
                    if not peptideList:
                        for sequence in sequence_list:
                            scores.append(predictor.predict_sequence(sequence,pred))
                    else:
                        scores.append(predictor.predict_sequence(sequence_list,pred,peptideList))
                    results.append((length, allele, scores, method_name))
            else:
                '''
                1. If the user supplied his own mhc sequence:
                2. Loop over sequences:
                    For each sequence, make predictions using 'netmhcpan'?
                '''
                mhc = 'User defined'
                length = (self.length).pop(0)
                predictor.initialize(mhc, int(length), self.hla_seq)
                method_name = "netmhcpan"
                scores = []
                if not peptideList:
                    for sequence in sequence_list:
                        scores.append(predictor.predict_sequence(sequence,pred))
                else:
                    scores.append(predictor.predict_sequence(sequence_list,pred,peptideList))
                results.append((int(length), mhc, scores, method_name))
        return results
     
    def method_lookup(self, allele, length):
        """TODO: extract all the references to the IEDB directories so that paths can be specified somewhere else"""
        import sqlite3 as lite
        path_to_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../iedb/mhc_i/data/allele.sqlite')
        conn =  lite.connect(path_to_db)
        method_list = []
        with conn:
            cur = conn.cursor()
            cur.execute("SELECT distinct(m.name) FROM mhci_allele a, mhci_method m \
            where a.method_id = m.id and a.name = ? and a.length = ? and \
            a.method_id != ? and a.method_id != ? and a.method_id != ? and a.method_id != ? and a.method_id != ? and a.method_id != ?", (allele, length, 1, 2, 5, 8, 9, 10))
            rows = cur.fetchall()
            methods = list(sum(rows,()))
            if len(methods) > 1 and 'netmhcpan' in methods:
                methods.remove('netmhcpan')
            method_list.extend(methods)
        return method_list
    
