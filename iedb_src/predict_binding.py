#!/usr/bin/env python

import sys
from optparse import OptionParser

from util import * #@UnusedWildImport
from seqpredictor import MHCBindingPredictions
from setupinfo import * #@UnusedWildImport

import pandas as pd

class Prediction():
    
    def __init__(self):
        self.row_data = []
        
    #TODO: needs to change 
    def read_protein(self, fname):
        f = open(fname, 'r')
        protein = Proteins()
        protein.extract_fasta(f.read())
        f.close()
        return protein
    
    def insert_dash(self, method_list, score_list):
        scores = []
        dash = "-","-"
        for score in score_list:
            score = list(score)
            if "smm" in method_list and not "ann" in method_list:
                score.insert(1,dash)
                score.insert(4,dash)
                score.insert(5,dash)
            if "ann" in method_list and "smm" in method_list and "comblib_sidney2008" in method_list:
                score.append(dash)
            if "comblib_sidney2008" in method_list and not "ann" in method_list and not "smm" in method_list:
                score.insert(1,dash)
                score.insert(2,dash)
                score.insert(5,dash)
            if "ann" in method_list and "smm" in method_list and not "comblib_sidney2008" in method_list:
                score.insert(5,dash)
                score.insert(6,dash)
            if "ann" in method_list and not "smm" in method_list and not "comblib_sidney2008" in method_list:
                score.insert(5,dash)
                score.insert(6,dash)
                score.insert(7,dash)
            if "netmhcpan" in method_list and not "smm" in method_list and not "ann" in method_list:
                score.insert(1,dash)
                score.insert(2,dash)
                score.insert(3,dash)
            scores.append(tuple(self.flatten(score)))
        return scores
    def predict(self, method,allele, length,peptides):
        df = self._predict(method,allele,length,peptides,peptideList = True)
        return df

    def _predict(self,method,allele,length,proteins,peptideList):
        species = get_species(allele)
        allele = allele.split()
        length = length.split()
        species = species.split()
        hla_seq = ''
        
        use_cutoff = cutoff_value = None

        input = InputData(self.version, method, allele, hla_seq, length, proteins, species)
        mhc_predictor = MHCBindingPredictions(input)
        if not peptideList:
            mhc_scores = mhc_predictor.predict(input.input_protein.sequences, peptideList = peptideList)
        else:
            mhc_scores = mhc_predictor.predict(proteins, peptideList = peptideList)
        
        table_rows = self.format_binding(input, mhc_scores, method, use_cutoff, cutoff_value, peptideList=peptideList, peptide_list=proteins)
        
        method_used = ','.join(mhc_predictor.get_method_set_selected(method))
        table_rows.sort(key=lambda tup: tup[6])
        table_rows = self.add_method_used(table_rows, method)
        
        # remove empty columns (netMHCpan related) from the result 
        if method == 'consensus':
            table_rows = [tuple(row for row in rows if row != '-') for rows in table_rows]
        
        # headers for different methods
        if method == 'IEDB_recommended':
            header = ('allele','seq_num','start','end','length','peptide','method',mhc_predictor.get_score_unit(),'ann_ic50','ann_rank','smm_ic50','smm_rank','comblib_sidney2008_score','comblib_sidney2008_rank','netmhcpan_ic50','netmhcpan_rank')
        elif method == 'consensus':
            header = ('allele','seq_num','start','end','length','peptide',mhc_predictor.get_score_unit(),'ann_ic50','ann_rank','smm_ic50','smm_rank','comblib_sidney2008_score','comblib_sidney2008_rank')
        else:
            header = ('allele','seq_num','start','end','length','peptide',mhc_predictor.get_score_unit(),'rank')

        d = {h:[] for h in header}
        for row in table_rows:
            for h,v in zip(header,row):
                d[h].append(v)
        L = len(d['allele'])
        d = {h:v for h,v in d.items() if len(v)==L}
        cols = [h for h in header if h in d]
        return pd.DataFrame(d)[cols]

    def commandline_input(self, args):
        '''Make predictions given user provided list of sequences. The input sequence is in fasta format.'''
        (method, allele, length, fname, peptideList) = args
        if not peptideList:
            proteins = self.read_protein(fname)
        else:
            with open(fname, 'r') as fh:
                proteins = [p.strip() for p in fh.readlines()]

        print self._predict(method,allele,length,proteins,peptideList).to_string()
        
    def modify(self, lst):
        return[tuple(self.flatten(x)) for x in lst]
        
    def flatten(self, tup): 
        from itertools import chain
        return list(chain(*(i if isinstance(i, tuple) else (i,) for i in tup)))
 
    def format_binding(self, proteins, results, method, cutoff, value, peptideList=False, peptide_list=None):
        if peptideList:
            proteins = peptide_list

        for length, allele, score, method_list in results:
            if method == 'consensus' or (method == 'IEDB_recommended' and method_list != 'netmhcpan'):
                score_list = []
                for s in score:
                    ranks_scores = reduce(lambda x, y: x + y, s[1])
                    scores = zip(s[0], zip(*ranks_scores))
                    scores = self.insert_dash(method_list, self.modify(scores))
                    score_list.append(scores)
                self.add_rows_binding(allele, length, proteins, score_list, method_list, cutoff, value, peptideList)
            elif method == 'IEDB_recommended' and method_list == 'netmhcpan':
                score_list = []
                for s in score:
                    scores = zip(s[0], zip(*s[1]))
                    scores = self.cons_netmhcpan(scores)
                    scores = self.insert_dash(method_list, self.modify(scores))
                    score_list.append(scores)
                self.add_rows_binding(allele, length, proteins, score_list, method_list, cutoff, value, peptideList)
            else:
                self.add_rows_binding(allele, length, proteins, score, method_list, cutoff, value, peptideList)
        return self.modify(self.row_data)
    
    
    def add_rows_binding(self, allele, length, proteins, score_list, method_list, cutoff, value, peptideList):
        if not peptideList:
            for (i,(sequence, predictions)) in enumerate(zip(proteins.input_protein.sequences, score_list)):
                for (k, prediction) in enumerate(predictions):
                    sequence_source = "%s" %(i+1)
                    sequence_start = "%s" %(k + 1)
                    sequence_end = "%s" %(k + int(length))
                    scanned_sequence = sequence[k : k + length]
                    self.row_data.append((allele, sequence_source, sequence_start, sequence_end, length, scanned_sequence, prediction, method_list.replace(",","-")))
        else:
            for (k,(sequence, prediction)) in enumerate(zip(proteins, score_list[0])):
                sequence_source = "%s" %(k+1)
                sequence_start = "1"
                sequence_end = "%s" %(int(length))
                scanned_sequence = sequence
                self.row_data.append((allele, sequence_source, sequence_start, sequence_end, length, scanned_sequence, prediction, method_list.replace(",","-")))
    
    def cons_netmhcpan(self, scores):
        score_list = []
        for score in scores:
            lis = list(score)
            del lis[-1]
            item2 = list(score[1])
            item2.append(score[0])
            lis.append(tuple(item2))
            score_list.append(tuple(lis))
        return score_list
    
    def add_method_used(self, table_rows, method):
        formated_data = []
        for row in table_rows:
            lis = list(row)
            if method == 'IEDB_recommended':
                if '-' not in lis[-1]:
                    lis.insert(6, lis[-1])
                else:
                    lis.insert(6, "Consensus ("+lis[-1].replace("-","/")+")")
                del lis[-1]
                formated_data.append(tuple(lis))
            else: 
                del lis[-1]
                formated_data.append(tuple(lis))
        return formated_data
    def commandline_input_mhc(self, parser, fname):
        """TODO: this function is not available programatically or with a peptide list"""
        '''This version takes a file containing an mhc sequence as input.
           Make predictions given user provided list of sequences. The input sequence is in fasta format.'''
        (options, args) = parser.parse_args()
        
        if len(args) == 3:
            (method, length, fname) = args
        else: (method, length) = args
        
        if method != 'netmhcpan': parser.error('Only netmhcpan has the option to take user-provided mhc sequence.')
        mhc = 'User defined'
        
        mhc_fh = open(options.filename_mhc, 'r')
        hla_seq = ''.join(mhc_fh.readlines())
        allele = [''.join(mhc_fh.readlines())]
        length = length.split()
        proteins = self.read_protein(fname)
        species = ['human']
        
        input_data = InputData(self.version, method, allele, hla_seq, length, proteins, species)
        
        length = length[len(length) - 1]
        
        mhc_predictor = MHCBindingPredictions(input_data)
        mhc_scores = mhc_predictor.predict(input_data.input_protein.sequences)
        (results_peptide_length, results_allele, scores, method_used) = mhc_scores[0]
    
        header = ['allele', 'length', 'peptide', mhc_predictor.get_score_unit()]
        print '\t'.join(header)
        
        for seq_index in range(len(scores)):
            seq = proteins.sequences[seq_index]
            peptide_list = get_peptides(seq, int(length))
            seq_scores = scores[seq_index]
            for (peptide, score) in zip(peptide_list, seq_scores):
                row = [mhc, length, peptide, score]
                print '\t'.join(map(str,row))


    def query_yes_no(self, question, default="yes"):
        """Ask a yes/no question via raw_input() and return their answer.
    
        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).
    
        The "answer" return value is one of "yes" or "no".
        """
        valid = {"yes":True,   "y":True,  "ye":True,
                 "no":False,     "n":False}
        if default == None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)


    def commandline_mhc(self, argv):
        '''Return all available MHC molecules against which predictions can be made.'''
        (method, mhc) = argv
        ms = MethodSet()
        method_index = ms.get_method_index(method)
        
        if method_index == 3:
            print "'arb' has been removed from the list of available methods."
            sys.exit(1)
        
        mhc_list = get_mhc_list(method_index)
    
        print "List of available (MHC,PeptideLength) for "+method
        if method == 'netmhcpan' or method == 'IEDB_recommended':
            answer = self.query_yes_no("The list is very long. Do you still want to print it?")
            if answer == False:
                exit
            else: 
                print "Species", "\t", "MHC", "\t", "PeptideLength"
                for (species,mhc, peptide_length) in mhc_list:
                    if mhc is not None:
                        print species, "\t", ""+mhc+"", "\t", peptide_length
        else: 
            print "Species", "\t", "MHC", "\t", "PeptideLength"
            for (species,mhc, peptide_length) in mhc_list:
                if mhc is not None:
                    print species, "\t", ""+mhc+"", "\t", peptide_length
                
    def commandline_method(self):
        '''Return all available prediction methods.'''
        ms = MethodSet()
        method_list = ms.get_method_list(version=self.version)
        print "MHC-I prediction methods:"
        print "-------------------------"
        for method in method_list:
            print method
        print
    
    def commandline_help(self):
        print " _______________________________________________________________________________________________________________________"
        print "|***********************************************************************************************************************|"
        print "| * List all available commands.                                                                                        |"
        print "| ./src/predict_binding.py                                                                                              |"
        print "|_______________________________________________________________________________________________________________________|"
        print "| * List all available MHC-I prediction methods.                                                                        |"
        print "| ./src/predict_binding.py method                                                                                       |"
        print "|_______________________________________________________________________________________________________________________|"
        print "| * List all available (MHC,peptide_length) for a given method.                                                         |"
        print "| ./src/predict_binding [method] mhc                                                                                    |"
        print "| Example: ./src/predict_binding.py ann mhc                                                                             |"
        print "|_______________________________________________________________________________________________________________________|"
        print "| * Make predictions given a file containing a list of sequences.                                                       |"
        print "| ./src/predict_binding [method] [mhc] [peptide_length] [input_file]                                                    |"
        print "| Example: ./src/predict_binding.py ann HLA-A*02:01 9 ./examples/input_sequence.fasta                                   |"
        print "|_______________________________________________________________________________________________________________________|"
        print "| * Make predictions given a file containing a list of sequences AND user-provided MHC sequence.                        |"
        print "| ** Only netmhcpan has this option.                                                                                    |"
        print "| ./src/predict_binding [method] -m [input_file_mhc] [peptide_length] [input_file]                                      |"
        print "| Example: ./src/predict_binding.py netmhcpan -m ./examples/protein_mhc_B0702.fasta 9 ./examples/input_sequence.fasta   |"
        print "|_______________________________________________________________________________________________________________________|"
        print "| * You may also redirect (pipe) the input file into the script.                                                        |"
        print "| Examples:                                                                                                             |"                                                                     
        print "| echo -e ./examples/input_sequence.fasta | ./src/predict_binding.py ann HLA-A*02:01 9                                  |"
        print "| echo -e ./examples/input_sequence.fasta | ./src/predict_binding.py netmhcpan -m ./examples/protein_mhc_B0702.fasta 9  |"
        print "|_______________________________________________________________________________________________________________________|"

    def main(self):
        import select
        self.version = '20130222'
        
        try:
            usage = "usage: %prog method allele or [options] arg length\n---\n\
Following are the available choices - \n\
   method: ann, comblib_sidney2008, consensus, IEDB_recommended, netmhcpan, smm, smmpmbec, pickpocket, netmhccons\n\
   allele: an allele name\n\
   length: a length"
            
            parser = OptionParser(usage=usage)
            parser.add_option("-m", dest="filename_mhc",
                              help="FILE containing a single MHC sequence in fasta format.", metavar="FILE")
            parser.add_option("-p", dest='peptide_list')
            (options, args) = parser.parse_args()
            
            # If there's input ready, do something, else do something
            # else. Note timeout is zero so select won't block at all.
            if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
                infile = sys.stdin.readline().strip()
                args.append(infile)
            
            if (len(args) == 0):                               self.commandline_help()
            elif ((len(args) == 1) and (args[0] == 'method')): self.commandline_method()
            elif ((len(args) == 2) and (args[1] == 'mhc')):    self.commandline_mhc(args)
            elif (len(args)  == 3):                            self.commandline_input_mhc(parser, args[2])  # args=[method, length, fname]
            elif (len(args)  == 4):                            self.commandline_input(args)  # args=[method, mhc, length, fname]
            else: 
                parser.error("incorrect number of arguments")
                self.commandline_help()
        except UnexpectedInputError,e:
            print str(e)
    def mymain(self):
        import argparse
        self.version = '20130222'

        parser = argparse.ArgumentParser(description='IEDB MHC-I binding prdiction tools.')

        parser.add_argument('--method', '-m', type=str,
                           choices = ['ann', 'comblib_sidney2008', 'consensus', 'IEDB_recommended', 'netmhcpan', 'smm','smmpmbec', 'pickpocket', 'netmhccons'],
                           help='method',
                           required=True)

        parser.add_argument('--allele','-a', type=str,
                           help='allele (use "list" to list available alleles for a given method)',
                           default='list')

        parser.add_argument('--filename_mhc', type=str,
                           help='file containing a single MHC sequence in fasta format')

        parser.add_argument('--length','-L',metavar='L', type=str, default = "9", help='length of the peptides')

        parser.add_argument('--fname', '-f', help='file containg amino acid sequences')

        parser.add_argument("--peptide_list","-p", default=False,help="input file contains a list of kmer peptides", action="store_true")

        args = parser.parse_args()

        if args.allele == 'list':
            self.commandline_mhc((args.method,''))
        elif not args.filename_mhc is None:
            self.commandline_input_mhc(args, (args.method,args.length,args.fname))  # args=[method, length, fname]
        else:
            self.commandline_input((args.method,args.allele,args.length,args.fname,args.peptide_list))  # args=[method, mhc, length, fname]



if __name__ == '__main__':
    Prediction().mymain()
    
