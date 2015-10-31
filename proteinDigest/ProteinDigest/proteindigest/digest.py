from Bio.SeqUtils.ProtParam import ProteinAnalysis

######################################################
######################################################
##                                                  ##
##  Base class is Digest                            ##
##  Specific enzyme classes inherit from Digest     ##
##      - specify cleavage sites and conditions     ##
##                                                  ##
##                                                  ##
######################################################
######################################################

#Base Class
class Digest(object):

    # factory is static method to create proper enzyme based on string passed
    def factory(type,peptide,misses,minlen,maxlen,minweight,maxweight):
        if type == "Chymotrypsin": return Chymotrypsin(peptide,misses,minlen,maxlen,minweight,maxweight)
        if type == "Chymotrypsin Low Specificity": return Chymotrypsin_low(peptide,misses,minlen,maxlen,minweight,maxweight)
        if type == "Trypsin": return Trypsin(peptide,misses,minlen,maxlen,minweight,maxweight)
        if type == 'Pepsin (ph >= 2.0)': return Pepsin_highph(peptide,misses,minlen,maxlen,minweight,maxweight)
        if type == 'Pepsin (ph = 1.3)': return Pepsin_lowph(peptide,misses,minlen,maxlen,minweight,maxweight)
        assert 0, "No enzyme " + type
    factory = staticmethod(factory)

    #Must give at least a peptide. Instantiates variables
    
    def __init__(self,peptide,misses=0,minlen=4,maxlen=1000,minweight=0,maxweight=999999999):
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
            if self.validCleave(i):
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
                        self.dpeps.append([l,len(l),ProteinAnalysis(str(l)).molecular_weight(),j,'',dp+r,str(1)+'-'+str(len(l))])
                if self.checkLenWeight(dp):
                    self.dpeps.append([dp,len(dp),ProteinAnalysis(str(dp)).molecular_weight(),j,l,r,str(self.sites[i]+2)+'-'+str(self.sites[i]+len(dp)+1)])
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

    # see interior comment for what each position of case corresponds to. Checks cases which make valid cleave false
    def validCleave(self,pos):
        #     0      1     2      3
        # (aaAtPos, pos, badAA, badAAPos)
        valid = False
        for case in self.cleaveSites:
            try:
                if self.peptide[pos+case[1]] in case[0]:
                    valid = True
                    if self.peptide[pos + case[3]] in case[2]:
                        return False
            except IndexError:
                continue
        return valid


##################################################################
##################################################################
##                                                              ##
##                      Enzyme Classes                          ##
##                                                              ##
##################################################################
##################################################################

## cleaveSites structure:
## (aa,positionOfAA,badAA,positionBadAA

class Chymotrypsin(Digest):
    cleaveSites = [
        ('W',0,'MP',1),
        ('YF',0,'P',1),
        ]


#Same as Chymotrypsin ^ but a few extra residues and criteria. Therefore just extend cleave sites.
class Chymotrypsin_low(Chymotrypsin):
    def __init__(self,peptide,misses=0,minlen=0,maxlen=1000,minweight=0,maxweight=999999999):
        self.cleaveSites.extend(
            [('M',0,'PY',1),
            ('H',0,'PDMW',1),
            ('L',0,'P',1)]
            )
        super(Chymotrypsin_low,self).__init__(peptide,misses,minlen,maxlen,minweight,maxweight)

class Trypsin(Digest):
    cleaveSites = [
        ('RK',0,'P',1),
        ]

class Pepsin_highph(Digest):
    cleaveSites = [
        ('FLWY',0,'RKH',-2),
        ('FLWY',0,'P',-1),
        ('FLWY',0,'P',2),
        ('FLWY',1,'RKH',-2),
        ('FLYW',1,'R',0),
        ('FLWY',1,'P',2),
        ('FLWY',1,'P',-1),
        ]

class Pepsin_lowph(Digest):
    cleaveSites = [
        ('FL',0,'RKH',-2),
        ('FL',0,'P',-1),
        ('FL',0,'P',2),
        ('FL',1,'RKH',-2),
        ('FL',1,'R',0),
        ('FL',1,'P',2),
        ('FL',0,'P',-1),
        ]
