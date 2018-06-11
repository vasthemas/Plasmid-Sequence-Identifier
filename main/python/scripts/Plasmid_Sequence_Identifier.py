#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 12:01:53 2018

@author: Vasanth Chandrasekhar
"""
import os
from itertools import chain
import xlwt
from Bio import SeqIO
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord


class plasmidIDtool(): 
    """
    Tool to search sequence of DNA and find regions of interest
    Support only .seq files
    """
    def __init__(self):
        
        self.wd = ''
        self.lis =[]
        self.region = ''
        self.searchRange = 0
        self.searchRegionLen = 0
        self.name =''
        self.complete = False 
        self.filename = ''
             
            
    def fetch_wd(self,wd): 
        """
        verify working dir exsists
        """
        try:
            os.chdir(wd)
        except OSError:
            print """Error: Sorry can't find this folder, recheck path"""
            return False
            
        else:
            print "Written content in the file successfully"
            self.wd = wd 
            return True
        
        
    def check_files(self):
        
        """
        generate list of .seq files 
        """
        self.lis = []
        direc = os.listdir(self.wd)
        
        for file_ in direc: 
            if file_.endswith(('.gb','.fa','.dna','.seq')): self.lis.append(file_)
                
        return self.lis 
    
    
    
    
    def ProduceCleanSeq(self,filename): 
        #MAKE SURE USELESS IS A CONSTANT 
        """
        Format file, removing unessential characters. Useless constant is the character length of 
        the file for .seq,.fa,.gb,.dna files.
    
        """
       
        if filename.endswith('.fa'):
            for seq_record in SeqIO.parse(filename, "fasta"): Data_formated = str(seq_record.seq).upper()
        elif filename.endswith('.gb'):
            for seq_record in SeqIO.parse(filename, "genbank"): Data_formated = str(seq_record.seq).upper()    
        elif filename.endswith('.dna'): 
            Data_formated = str(snapgene_file_to_dict(filename)['seq']).upper()
        elif filename.endswith('.seq'):
                useless = 12 
                data = open(filename).read().split()
                Data_unformatted = list(chain.from_iterable(data[i] for i in range(len(data))))
                Data_formated = Data_unformatted[useless:]
                
        Data_str = ''.join(Data_formated)
        Clean_Seq = ''.join([i for i in Data_str if not i.isdigit()])
        
        return Clean_Seq


    def find_region(self,Region,clean_seq):
        """
        return bool if region is in sequence
        """
       
        return  Region in clean_seq
    
    
    
    def find_region_range(self,Region,clean_seq,searchRange,searchRegionLen): 
        #DEAL WITH THE WRAPAROUND ISSUE 
        
        """
        Look for given region in the cleaned sequence
        Foreward is 3' to 5'
        """
        
        
        if self.find_region(Region,clean_seq) == False: 
            print "Region does not exsist in given sequence "
            print "  "
            new_regF, new_regB = "NaN","NaN"
        
        else:
            print "Region Found"
            print "  "
            regionlen = len(Region)
            region_start = clean_seq.find(Region)
            region_end = region_start + regionlen
            new_pos_F = region_end + searchRange
            new_pos_B = region_start - searchRange
            
            #deal with wrap around issue. 
            if new_pos_B < 0: 
                #THINK OF BETTER WAY TO DO THIS
                new_pos_B = len(clean_seq) * 2
            
            new_regF = clean_seq[new_pos_F : new_pos_F + searchRegionLen]
            new_regB = clean_seq[new_pos_B - searchRegionLen: new_pos_B]

        return new_regF, new_regB

   

    
    def rev_comp(self,seq): 
        """
        Produce the reverse complement sequence of given sequence
        """
        if seq == "NaN": rev_comp = "NaN"
        else:     
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rev_comp = ''.join([complement[base] for base in seq[::-1]])
        
        return rev_comp
    
    def DNA_to_protein(self,seq): 
        """
        Convert DNA to Amino acid
        """
        if seq == "NaN": protein = "NaN" 
        elif "N" in seq: protein = "Invalid sequence"
       
        else:
            table = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
            }

            protein =""
            if len(seq)%3 == 0:
                for i in range(0, len(seq), 3):
                    codon = seq[i:i + 3]
                    protein+= table[codon]
        return protein 

    
    def start(self): 
        
        """
        Runs tool UI
        """
        while self.complete == False:
            
            #Enter the path name by user 
            wd = str(raw_input('enter the folder path you would like to analyze: '))
            
            if self.fetch_wd(wd) == True: 
                print 'Checking files'
                self.check_files()
                print "You have" + " " + str(len(self.lis)) + " " + "file(s)"
                for  i in self.lis:
                    print i 
            else: 
                print "  "
                print "Cant Find you File, lets try again"
            
            if len(self.lis) == 0: print "nothing in here, lets try again"    
            else:
                #Enter specifics for user
                #self.region = str(raw_input('enter the refernce region you would like to search: '))
                self.region = (str(raw_input('enter the refernce region(s) you would like to search, seperate regions with space: '))).split()
                self.searchRange = int(raw_input('how many base pairs would you like to search : '))
                self.searchRegionLen = int(raw_input('how many base pairs is region you like to find : '))
                self.name = raw_input('what you like to call your excel file : ')
                
                #Generate Excel files 
                 #DONT NEED BACKWARD SEARCHER 
                row = 0
                book = xlwt.Workbook(encoding="utf-8")
                sheet1 = book.add_sheet("Data")
                sheet1.write(0, 0, "File Name")
                sheet1.write(0, 1, "Reference Region")
                sheet1.write(0, 2, "Search Distance")
                sheet1.write(0, 3, "New Region Foreward")   
                sheet1.write(0, 4, "New Region Backward")
                sheet1.write(0, 5, "Protein Foreward")
                sheet1.write(0, 6, "Protein Backward")
                sheet1.write(0, 7, "New Region Foreward(Rev Comp)")
                sheet1.write(0, 8, "New Region Backward(Rev Comp)")
                sheet1.write(0, 9, "Protein Foreward(Rev Comp)")
                sheet1.write(0, 10, "Protein Backward(Rev Comp)")
                
                for reg in self.region:
                    for idx in range(len(self.lis)):
                        print self.lis[idx]
                        #Correct for row in excel
                        row += 1 
                       
                        #Clean up sequence 
                        clean_seq = self.ProduceCleanSeq(self.lis[idx])

                        #Find regions and convert to proteins 
                        new_regF, new_regB = self.find_region_range(reg,clean_seq,self.searchRange,self.searchRegionLen)
                        proteinF = self.DNA_to_protein(new_regF)
                        proteinB = self.DNA_to_protein(new_regB)
                        new_regF_rc = self.rev_comp(new_regF)
                        new_regB_rc = self.rev_comp(new_regB)
                        proteinF_rc = self.DNA_to_protein(new_regF_rc)
                        proteinB_rc = self.DNA_to_protein(new_regB_rc)

                        #enter in the file 
                        sheet1.write(row, 0, self.lis[idx])
                        sheet1.write(row, 1, reg)
                        sheet1.write(row, 2, self.searchRange)
                        sheet1.write(row, 3, new_regF)
                        sheet1.write(row, 4, new_regB)
                        sheet1.write(row, 5, proteinF)
                        sheet1.write(row, 6, proteinB)
                        sheet1.write(row, 7, new_regF_rc)
                        sheet1.write(row, 8, new_regB_rc)
                        sheet1.write(row, 9, proteinF_rc)
                        sheet1.write(row, 10,proteinB_rc)

                self.complete = True 
                book.save(self.name)
                print "file saved"
                return self.complete
                
                

                
                
                
def main():
    plasmidIDer = plasmidIDtool()
    plasmidIDer.start()

if __name__ == '__main__':
    main()
 
