

import os
from itertools import chain
import xlwt
from Bio import SeqIO
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord

class SequenceGrabber(): 
    """
    Tool to search sequence of DNA to find insert between given primer regions,
    Translates insert to protein.
    Converts insert to reverse complement.
    Compiles all data into an excel file.
    Supports only .seq, .gbk, .fa, .dna
    """
    def __init__(self):
        
        self.wd = ''
        self.lis =[]
        self.regionup = ''
        self.regiondown = ''
        self.searchRange = 0
        self.searchRegionLen = 0
        self.name =''
        self.complete = False 
        self.filename = ''
             
            
    def fetch_wd(self,wd): 
        """
        verify working dir exists
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
        generate list of supported files
        """
        self.lis = []
        direc = os.listdir(self.wd)
        
        for file_ in direc: 
            if file_.endswith(('.gb','.fa','.dna','.seq')): self.lis.append(file_)
                
        return self.lis 
    
    
    def ProduceCleanSeq(self,filename):
        """
        Format file, removing unessential characters.
    
        """
        if filename.endswith('.fa'):
            for seq_record in SeqIO.parse(filename, "fasta"): Data_formatted = str(seq_record.seq).upper()
        elif filename.endswith('.gb'):
            for seq_record in SeqIO.parse(filename, "genbank"): Data_formatted = str(seq_record.seq).upper()    
        elif filename.endswith('.dna'): 
            Data_formatted = str(snapgene_file_to_dict(filename)['seq']).upper()
        elif filename.endswith('.seq'):
                
                data = open(filename).read().split()
                data.pop(0)
                Data_formatted = list(chain.from_iterable(data[i] for i in range(len(data))))

        Data_str = ''.join(Data_formatted)
        Clean_Seq = ''.join([i for i in Data_str if not i.isdigit()])
        
        return Clean_Seq


    def find_region(self,Region,clean_seq):
        """
        return bool if region is in sequence
        """
        return  Region in clean_seq
    
    
    def grab_sequence(self,Regionup,Regiondown,clean_seq):
        """
        Find the insert between the given upstream and downstream primers
        on the formatted sequence.

        """

        if (self.find_region(Regionup,clean_seq) == False and self.find_region(Regiondown,clean_seq) == False): 
            print "Region does not exsist in given sequence "
            print ' '
            seq = "NaN"
            return seq
            
        else: 
            print "Regions Found"
            print "  "
            regiondown_start = clean_seq.find(Regiondown)
            regionup_start = clean_seq.find(Regionup) 
            regionupL = len(Regionup)
            regiondownL = len(Regiondown)
            
            #confirm region order, reverse if incorrect 
            if  regiondown_start <  regionup_start:
                regiondown_start = clean_seq.find(Regionup)
                regionup_start = clean_seq.find(Regiondown)
                regionupL = len(Regiondown)
                regiondownL = len(Regionup)
                
                print 'given down stream region is upstream of given upstream region'
                print ' ' 
    
            regionup_end = regionup_start + regionupL
            
            seq = clean_seq[regionup_end: regiondown_start]
            
            return seq
            
    def find_region_range(self,Region,clean_seq,searchRange,searchRegionLen):
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
        if 'N' in seq: rev_comp = "NaN"
        else:     
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            rev_comp = ''.join([complement[base] for base in seq[::-1]])
        
        return rev_comp
    
    
    def DNA_to_protein(self,seq): 
        """
        Convert DNA to Amino acid, translation continues until out of frame
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

            off_frame = len(seq)%3
            for i in range(0, len(seq) - off_frame, 3):
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
                for file_ in self.lis:
                    print file_
                    
            else: 
                print "  "
                print "Cant find your file, lets try again"
            
            if len(self.lis) == 0: print "nothing in here, lets try again"    
            
            else:
                #Enter specifics for user

                self.regionup,self.regiondown = (str(raw_input('enter the refernce region(s) you would like to search(upstream,downstream), seperate regions with space: '))).split()
                self.name = raw_input('what you like to call your excel file : ')
                
                #Generate Excel files 
                row = 0
                book = xlwt.Workbook(encoding="utf-8")
                sheet1 = book.add_sheet("Data")
                sheet1.write(0, 0, "File Name")
                sheet1.write(0, 1, "Reference Region upstream")
                sheet1.write(0, 2, "Reference Region downstream")
                sheet1.write(0, 3, "Insert")
                sheet1.write(0, 4, "Insert (reverse complement)")
                sheet1.write(0, 5, "Translated Insert")
                sheet1.write(0, 6, "Translated Insert (reverse complement)")
                
                
#                 for reg in self.region:
                for idx in range(len(self.lis)):
                    print self.lis[idx]
                    #Correct for row in excel
                    row += 1 
                    
                    #Clean up sequence 
                    clean_seq = self.ProduceCleanSeq(self.lis[idx])

                    seq = self.grab_sequence(self.regionup,self.regiondown,clean_seq)
                    rev_comp_seq = self.rev_comp(seq)

                    protein = self.DNA_to_protein(seq)
                    protein_rec_comp = self.DNA_to_protein(rev_comp_seq)


                    sheet1.write(row, 0, self.lis[idx])
                    sheet1.write(row, 1, self.regionup)
                    sheet1.write(row, 2, self.regiondown)
                    sheet1.write(row, 3, seq)
                    print 'Insert'
                    print seq
                    print " "
                    sheet1.write(row, 4, rev_comp_seq)
                    sheet1.write(row, 5, protein)
                    sheet1.write(row, 6, protein_rec_comp)
                    
                
                    

                self.complete = True 
                book.save(self.name)
                print "file saved"
                return self.complete
                
                

                
                
                
def main():
    plasmidIDer = SequenceGrabber()
    plasmidIDer.start()

if __name__ == '__main__':
    main()
 
            

