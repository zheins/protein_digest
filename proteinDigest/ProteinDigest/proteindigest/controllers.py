# symbols which are imported by "from proteindigest.controllers import *"
__all__ = ['Root']

from turbogears import controllers, expose, flash
from turbogears import validate, validators
from turbogears import widgets, error_handler
from digest import *
from Bio.Seq import Seq
import Bio.SeqIO
import sys,Bio.Alphabet,Bio.Alphabet.IUPAC as IUPAC
import urllib
import re


# Custom validator, validates whether Uniprot accession or proper aa seq
class accesionOrSeq (validators.FancyValidator):
    def _to_python(self, value, state):
        aos = value.strip()

        badAccession = False
        badSeq = False
        fasta = False

        #First check if accession. Check by attempting to get fasta from Uniprot. If fails, badAccession will be set to true
        url = "http://www.uniprot.org/uniprot/" + value.strip() + ".fasta"
        handle = urllib.urlopen(url)
        if aos.startswith('>'):
            seq = ''.join(aos.splitlines()[1:])
            badAccession = True
            fasta = True
        else:
            seq = aos
            try:
                record = Bio.SeqIO.read(handle,'fasta')
            except ValueError:
                badAccession = True

        # If it's a bad accession, try it out as a then simple sequence
        if badAccession:
            seq = seq.upper()
            seq = re.sub(r'[^\x20-\x7E]','',seq)
            
    
            # Using IUPAC amino acid alphabet. Verify sequence only has these characters
            peptide = Seq(seq,IUPAC.protein)
            if not Bio.Alphabet._verify_alphabet(peptide):
                badSeq = True

        #If badAccession and badSeq, raise as invalid. Both must be bad.
        if badAccession and badSeq:
            raise validators.Invalid(
                'Invalid Uniprot Accession or Sequence',
                value, state)
        return value

# Using turbogears to set up input fields
class InputFields(widgets.WidgetsList):
    peptide = widgets.TextArea(label="Uniprot Accession or Seq")
    
    enzyme = widgets.SingleSelectField(label="Enzyme",
                                     options=["Chymotrypsin",
                                              "Chymotrypsin Low Specificity",
                                              "Trypsin",
                                              "Pepsin (ph = 1.3)",
                                              "Pepsin (ph >= 2.0)"],
                                     default="Trypsin")
    misses = widgets.TextField(label="Allowed Missed Cleavages")
    minlen = widgets.TextField(label="Minimum Peptide Length")
    maxlen = widgets.TextField(label="Maximum Peptide Length")
    minweight = widgets.TextField(label="Minimum Peptide Weight")
    maxweight = widgets.TextField(label="Maximum Peptide Weight")

    
       
# Turbogears, set validation schema
class InputFieldsSchema(validators.Schema):

    # Regex validator ensures only certain characters are input, essentially alphanumeric with new lines and spaces
    peptide = validators.All(accesionOrSeq(),validators.String(min=3),validators.Regex(r'^[^\x00-\x09\x0B-\x0C\x0E-\x1F\x21-\x2F\x3A-\x3C\x3F-\x40\x5B-\x5E\x60\x7B\x7D-\x7F]*$'))
    enzyme = validators.OneOf(["Chymotrypsin","Chymotrypsin Low Specificity","Trypsin","Pepsin (ph = 1.3)","Pepsin (ph >= 2.0)"])
    misses = validators.Int(if_empty=0,min=0, max=50)
    minlen = validators.Int(if_empty=0,min=0)
    maxlen = validators.Int(if_empty=1000000000,min=0)
    minweight = validators.Int(if_empty=500,min=0)
    maxweight = validators.Int(if_empty=1000000000,min=0)
    
# Turbogears form creation using input fields and validation schema
input_form = widgets.TableForm(
    fields = InputFields(),
    action = "digest_submit",
    submit_text = "Digest",
    validator = InputFieldsSchema()
    )
    
class Root(controllers.RootController):
    """The root controller of the application."""

    @expose(template="proteindigest.templates.welcome")
    def index(self):
        """"Show the welcome page."""
        # log.debug("Happy TurboGears Controller Responding For Duty")
        return dict(form=input_form)

    
    @expose(template="proteindigest.templates.digest_submit")
    @expose(template="proteindigest.templates.digest_submitxml",as_format='xml',format='xml')
    @validate(form=input_form)
    @error_handler(index)
    def digest_submit(self,peptide,enzyme,misses,minlen,maxlen,minweight,maxweight):
        accession = True
        # Try to get fasta from uniprotr
        url = "http://www.uniprot.org/uniprot/" + peptide.strip() + ".fasta"
        handle = urllib.urlopen(url)
        if peptide.startswith('>'):
                peptide = ''.join(peptide.splitlines()[1:])
        else:
            try:
                record = Bio.SeqIO.read(handle,'fasta')
                peptide = record.seq
                pepname = record.name

            # If this block executes, was not a valid uniprot entry, handle as pure sequence
            except ValueError:
                accession = False
        if not accession:       
            peptide = peptide.upper()
            peptide = re.sub(r'[^\x41-\x5A]','',peptide)
            peptide = Seq(peptide,IUPAC.protein)
        
        # Instantiate Digest object using factory, passing enzyme for type and all inputs    
        digester = Digest.factory(enzyme,peptide,misses,minlen,maxlen,minweight,maxweight)

        # Use the digest() method to get the peptides
        peptides = digester.digest()

        # This is purely formatting for output, inserting commas between sites
        formattedSites = ', '.join(map(str,map(lambda x: x+1,digester.sites)))

        # Chops off last comma
        formattedSites = formattedSites[:len(formattedSites)]

        lines = []

        # print out numbers on top of corresponding char
        # Printing 60 amino acids per line, with position at begining and end ( 1 MGC... 60)
        # a string of spaces and asterisks is created according to what sites have cleaves, placed above the amino acid string
        for i in range(0,len(peptide),60):
            if (i+60) < len(peptide):
                lines.append('\t' + ''.join(map(lambda x: '*' if x in digester.sites else ' ',range(i,i+60))) \
                             + '\n' + str(i+1) + '\t' + peptide[i:i+60] + '\t' + str(i+60))
            else:
                lines.append('\t' + ''.join(map(lambda x: '*' if x in digester.sites else ' ',range(i,i+60))) + '\n' + str(i+1) + '\t' + peptide[i:i+60] + '\t' + str(len(peptide)))
    
        return dict(peps=peptides,di=digester,enzyme=enzyme,lines=lines,fs=formattedSites)
