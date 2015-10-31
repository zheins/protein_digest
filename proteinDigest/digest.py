from Bio.SeqUtils.ProtParam import ProteinAnalysis

######################################################
######################################################
##                                                  ##
##  Base class is Digest                            ##
##  Specific enzyme classes inherit from Digest     ##
##      - specify cleavage sites                    ##
##      - define criteria for cleavage              ##
##                                                  ##
######################################################
######################################################

#Base Class
class Digest(object):

    #Must give at least a peptide. Instantiates variables
    
    def __init__(self,peptide,misses=0,minlen=0,maxlen=1000,minweight=0,maxweight=999999999):
        self.peptide = peptide
        self.misses = int(misses)
        self.minlen = int(minlen)
        self.maxlen = int(maxlen)
        self.minweight = int(minweight)
        self.maxweight = int(maxweight)
        self.sites = []

        #Will be a list of lists, each list will be a peptide with the following structure:
        #Peptide, length, weight, missed cleavages, amino acids to left, amino acids to right
        
        self.dpeps = []

    #Cleaves, analyzes, and returns list of peptides (see description in __init__)
        
    def digest(self):
        self.cleave()
        self.analyzeCleaves()
        return self.dpeps

    #Get sites of all possible cleavage locations in peptide
    
    def cleave(self):
        sites = []

        #For each residue, if it a cleavable amino acide and meets all criteria, add it to the list of cleavage sites
        
        for i,residue in enumerate(self.peptide):
            if residue in self.cleaveSites and self.validCleave(i):
                sites.append(i)
        self.sites = sites

    #For each peptide specified by cleave(), determines if it meets user specified criteria and fills dpeps as specified in __init__
        
    def analyzeCleaves(self):

        #i used to iterate through cleave sites
        #j used to iterate for miss cleaves. Skips j cleave site(s) when calculating the peptide from cleave sites
        
        for i in range(len(self.sites)):
            end = False
            for j in range(self.misses+1):
                l = self.peptide[:self.sites[i]+1]
                try:
                    r = self.peptide[self.sites[i+j+1]+1:]
                    dp = self.peptide[self.sites[i]+1:self.sites[i+j+1]+1]
                except IndexError:
                    #When code reaches this block, it means the end of the input string has been found
                    #Set end to true to stop going through missed cleaves, no more exist
                    r = ''
                    dp = self.peptide[self.sites[i]+1:]
                    end = True
                if i == 0:
                    l = self.peptide[:self.sites[i+j]+1]
                    if self.checkLenWeight(l):
                        self.dpeps.append([l,len(l),ProteinAnalysis(str(l)).molecular_weight(),j,'',dp+r])
                if self.checkLenWeight(dp):
                    self.dpeps.append([dp,len(dp),ProteinAnalysis(str(dp)).molecular_weight(),j,l,r])
                if end:
                    break
            

    #returns True or False based on length and weight of resulting peptide
                
    def checkLenWeight(self,pep):
        if len(pep) < self.minlen or len(pep) > self.maxlen:
            return False
        else:
            weight = ProteinAnalysis(str(pep)).molecular_weight()
            if weight < self.minweight or weight > self.maxweight:
                return False
        return True



##################################################################
##################################################################
##                                                              ##
##                      Enzyme Classes                          ##
##                                                              ##
##################################################################
##################################################################
        
class Chymotrypsin(Digest):
    badP1prime = 'MP'
    cleaveSites = ['W','Y','F']
    def validCleave(self,pos):
        if self.peptide[pos+1] in self.badP1prime:
            return False
        else:
            return True

#Same as Chymotrypsin ^ but a few extra residues and criteria
class Chymotrypsin_low(Chymotrypsin):
    badP1primeM = 'Y'
    badP1primeH = 'NMW'
    def __init__(self,peptide,misses=0,minlen=0,maxlen=1000,minweight=0,maxweight=999999999):
        self.cleaveSites.extend(['M','L','H'])
        super(Chymotrypsin_low,self).__init__(peptide,misses,minlen,maxlen,minweight,maxweight)
        
    def validCleave(self,pos):
        #use Chymotrypsin validCleave to check its validity, then check for low specificity validity
        vc = super(Chymotrypsin_low,self).validCleave(pos)
        if (not vc) or (self.peptide[pos] == 'M' and self.peptide[pos+1] in self.badP1primeM) or (self.peptide[pos] == 'H' and self.peptide[pos+1] in self.badP1primeH):
            return False
        else:
            return True

class Trypsin(Digest):
    badP1prime = 'P'
    cleaveSites = ['R','K']

    def validCleave(self,pos):
        if self.peptide[pos+1] in self.badP1prime and not ((self.peptide[pos] == 'K' and self.peptide[pos-1] == 'W') or (self.peptide[pos] == 'R' and self.peptide[pos-1] == 'M')):
            return False
        else:
            return True
    
    
