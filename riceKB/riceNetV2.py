import sys
from riceKB.globalVars import *
from riceKB.utils import *
import pprint
import re
import os
import pandas as pd
import numpy as np

'''
Created on Dec, 2019
The interpro module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling interpro data
It runs with ParentChildTreeFile.txt file downloaded from EBI interpro website

1 - set up the file input argument

2 - run the program 

@author: larmande
'''
__author__  = "larmande"


def goldStandard(infile,outfile):
    #array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str')
    with open(infile, 'r', encoding='utf-8') as file:
        rdf_writer = open(outfile, "w")
        nb_line = 0
        rdf_writer.write(str(getRDFHeaders()))
        for line in file:
            buffer = ''
            gene1, gene2 = line.split('\t')
            gene2 = re.sub('\s','',gene2)
            buffer += "<" + base_resource_uri  + gene1 + ">"
            buffer += "\t" + obo_ns + "RO_0002434" + "\t"
            buffer += "<" + base_resource_uri  + gene2 + "> ;\n"
            buffer = re.sub(' ;$', ' .\n', buffer)
            rdf_writer.write(buffer)
            print(buffer)
        print("nb de lignes traitees:" + str(nb_line))

path = '/Users/plarmande/workspace2015/datasets/RiceNet_GoldSandard.txt'
path_output = '/Users/plarmande/workspace2015/datasets/RiceNet_GoldSandard.ttl' # The output

#pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

#RDFConverter(ds, path_output)

goldStandard(path,path_output)