from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
import sys
#print sys.path
from riceKB.globalVars import *
import pprint
import re
import os
import pandas as pd
import numpy as np
import urllib
import json

'''
Created on Sept, 2019
The ensmbl.rdf module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Esembl RDF datasets to normalize with AgroLD

TODO:
    1) Add documentation
    2) better Error handling
    
Questions: 
- which URI to define for resources such as genes, transcript and proteins ? \
for now its all after http://www.southgreen.fr/agrold/resource/ 

@author: larmande
'''

#todo modify term:Gene to agrold_vocab:Gene > check all the terms use in Ensembl vocabulary and aling with agrold vocabulary
#todo check the chromosome URI and decide how to modify

__author__  = "larmande"

prefixes = {}
patern1 = "^<http://rdf.ebi.ac.uk/resource/ensembl/>"
def ensemblParser(files):
    """
    parse the ensembl rdf file
    :param files: path of the file
    :return: python objects corresponding to prefix and triples
    """
    if(os.path.isfile(str(files))):
        print("***************** Parsing Esembl RDF data ********************\n")
        with open(files) as fp:
            for line in fp:
                if(re.match(r'^@prefix', line)):
                    #print line
                    prefix, value = re.split(':\s+', line)
                    value = re.sub('\.\n$', '', value)
                    prefixes[prefix]=value
                else:
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl/\d*/?', 'http://www.southgreen.fr/agrold/resource/', line)
                    if re.findall("term:inEnsemblAssembly",line):
                        line = re.sub('term:inEnsemblAssembly',
                                      'term:inAssembly', line)
                    if re.findall("term:inEnsemblSchemaNumber",line):
                        line = re.sub('term:inEnsemblSchemaNumber',
                                      'term:inSchemaNumber', line) # a supprimer ou modifier
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.transcript/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.transcript/',
                                      '<http://www.southgreen.fr/agrold/resource/transcript/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.protein/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.protein/',
                                      '<http://www.southreen.fr/agrold/resource/protein/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.exon/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.exon/',
                                      '<http://www.southreen.fr/agrold/resource/exon/', line)
                    if re.findall("term:protein",line):
                        line = re.sub('term:protein', 'term:Protein', line)
                    if re.findall("term:Protein_coding",line):
                        line = re.sub('term:Protein_coding', 'term:Gene', line)
                    if re.findall("term:protein_coding",line):
                        line = re.sub('term:protein_coding', 'term:Gene', line)
                    if re.findall("term:EnsemblRegion", line):
                        line = re.sub('term:EnsemblRegion', 'term:Region', line)
                    if re.findall("<http://rdf.ebi.ac.uk/dataset/ensemblgenomes/", line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/dataset/ensemblgenomes/\d*/?',
                                      '<http://www.southgreen.fr/agrold/dataset/ensemblgenomes/', line)
                    # change   obo:SO_has_part >  term:has_part
                    # obo:SO_transcribed_from > term:transcribed_from ou term:develops_from
                    # obo:SO_translates_to > term:translates_to ou term:encodes

                    print(line)




    else:
        print("***************** File not found ********************\n")
        print(files)


ROOT_DIR = '/Users/plarmande/Downloads/larmande/'

ensembl_files =  os.path.join(ROOT_DIR+'sample.os.ttl')
ensembl_out = '/Users/plarmande/Downloads/larmande/'


pp = pprint.PrettyPrinter(indent=4)

ensemblParser(ensembl_files)
print("***************** Esembl RDF data ********************\n")

print("********************************************************\n\n")

# if __name__ == '__main__':
# 	link = "http://data.gramene.org/v60/genes?q=Os08g0494375"
# 	result = connectionError(link)
# 	print(result.content)