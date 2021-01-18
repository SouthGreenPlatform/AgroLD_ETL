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
The interpro module is created as part of the AgroLD project.

This module contains Parsers, RDF converters and generic functions for handling interpro data
It runs with ParentChildTreeFile.txt file downloaded from EBI interpro website

1 - set up the file input arguments "path" for the input file and "path_output" for the output file

2 - run the program 
fileParser(path,path_output)

@author: larmande
'''
__author__  = "larmande"

def fileParser(infile,outfile):
    #array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str')
    with open(infile, 'r', encoding='utf-8') as file:
        root_item = ""
        current_item = ""
        second_item  = ""
        third_item = ""
        fourth_item = ""
        rdf_writer = open(outfile, "w")
        nb_line = 0
        rdf_writer.write(str(getRDFHeaders()))
        for line in file:
            buffer = ''
            if not line.startswith('--'):
                nb_line += 1
                key,value,undef = line.split('::')
                root_item = key
                buffer += "<" + interpro_ns  + key + ">\n"
                buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein_Domain" + " ;\n"
                buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + value + "\" ;\n"
                buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + key + "\" ;\n"
            elif not line.startswith('----'):
                nb_line += 1
                key, value, undef = line.split('::')
                key = re.sub('--','',key)
                second_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein_Domain" + " ;\n"
                buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + value + "\" ;\n"
                buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + key + "\" ;\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + "\"" + root_item + "\" ;\n"
            elif not line.startswith('------'):
                nb_line += 1
                key, value, undef = line.split('::')
                key = re.sub('--', '', key)
                third_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein_Domain" + " ;\n"
                buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + value + "\" ;\n"
                buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + key + "\" ;\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + "\"" + second_item + "\" ;\n"
            elif not line.startswith('--------'):
                nb_line += 1
                key, value, undef = line.split('::')
                key = re.sub('--', '', key)
                fourth_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein_Domain" + " ;\n"
                buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + value + "\" ;\n"
                buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + key + "\" ;\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + "\"" + third_item + "\" ;\n"
            else:
                nb_line += 1
                key, value, undef = line.split('::')
                key = re.sub('--', '', key)
                parent_item = current_item
                current_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein_Domain" + " ;\n"
                buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + value + "\" ;\n"
                buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + key + "\" ;\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + "\"" + fourth_item + "\" ;\n"

            buffer = re.sub(' ;$', ' .\n', buffer)
            rdf_writer.write(buffer)
            print(buffer)
        print("nb de lignes traitees:" + str(nb_line))



#TEST PARAM
path = '/Users/plarmande/workspace2015/datasets/ParentChildTreeFile.txt'
path_output = '/Users/plarmande/workspace2015/datasets/ParentChildTreeFile.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
#ds = geneParser(path)   # The parsing file
#pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

#RDFConverter(ds, path_output)

fileParser(path,path_output)