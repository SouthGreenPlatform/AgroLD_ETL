import sys
print (sys.path)
from riceKB.globalVars import *
from riceKB.globalVars import base_vocab_ns
from riceKB.gffParser import *
from riceKB.utils import *
import pprint
import re
import os
import json
import pandas as pd
import numpy as np
import rdflib
from rdflib.graph import Graph

# TODO modifier dans la base les predicats develops_from et has_trait avec developsFrom hasTrait
'''
Created on Jan, 2021
This module is created as part of the AgroLD Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Pearl Millet genome annotation \
datasets to normalize with AgroLD

TODO:
    1) Add documentation
    2) better Error handling

Questions: 
- which URI to define for resources such as genes, transcript and proteins ? \
for now its all after http://www.southgreen.fr/agrold/resource/ 

@author: larmande
'''


__author__  = "larmande"

ROOT_DIR = ''
file_input = ''
prefixes = {}

def milletParser(files, type):
    """
       parse the ensembl rdf file
       :param files: path of the file
       write  objects corresponding to prefix and triples
    """
    if (type == "annotation"):
        print('OK annotation')

        file1 = "annotation_modified.tsv"
        file_path = os.path.join(data_dir, file1)
        annotation_hash = {}

        with open(file_path, "r") as filehandle:
            for line in filehandle:
                line = line.rstrip('\n')
                line = re.sub('\"', '', line)

                records = line.split('\t')
                geneid = records.pop(0)
                db = records.pop(0)
                if geneid in annotation_hash:
                    annotation_hash[geneid].append({'type': db, 'value': records})
                else:
                    annotation_hash[geneid]=[]
                    annotation_hash[geneid].append({'type':db, 'value': records})
        return annotation_hash
    else:

        file1 = "liste_gene_autre.txt"
        path = os.path.join(data_dir, file1)
    #output_writer = open(output_file, "w")
        print("*************Parsing Pearl Millet gene list data ***********\n")
        array = pd.read_csv(path, sep='\t', delimiter=None, dtype='str',\
                                              skip_blank_lines=True)
        array.fillna('', inplace=True)
        print(array.head())
        row_count = array.shape[0]
        print("Number of genes: %s\n" % (str(row_count)))
        print(" gene list data has been parsed!\n")
        print("*************************************\n\n")
        return array

def annotation2RDF(annotation, output):
    print("************* RDF conversion begins***********\n")

    ttl_handle = open(path_output, "w")
    pub_handle = open(uniprotid_list, 'w')

    ttl_handle.write(str(getRDFHeaders()))

    for geneid in annotation:
        annotation_list = annotation[geneid]
        #print(len(annotation_list))
        for i in range(len(annotation_list)):
            annotation_dict = annotation[geneid][i]
            print(annotation_dict['type'])
            if (annotation_dict['type']== 'GO'):
                print(annotation_dict['value'])
            elif (annotation_dict['type']== 'InterPro'):
                print(annotation_dict['value'])
            elif (annotation_dict['type']== 'KEGG'):
                print(annotation_dict['value'])
            elif (annotation_dict['type']== 'Swiss-Prot'): # InterPro #KEGG
                print(annotation_dict['value'])
            elif (annotation_dict['type'] == 'TrEMBL'):
                print(annotation_dict['value'])
            #for key in annotation_dict.keys():
                #print(annotation_dict['value'])
            # for annotation_dict in annotation[geneid][i]:
            #     print(annotation_dict)
            #     if annotation_dict['type']=='GO':
            #         print(annotation_dict['value'])
    # for key in annotation.keys():
    #     taxon_id = "4543"
    #     print(key)
    #     for value_type in annotation[key]['type']:
    #         if annotation[key]['type'] != 0:
    #             print(value_type)


#TEST PARAM

#path = '/Users/pierre/workspace2015/datasets/Pearl_Millet/liste_gene_autre.txt'
rdf_file = 'annotation.ttl' # The output
#pub_dict = '/Users/plarmande/workspace2015/datasets/uniprot-dictionnary-millet.txt'
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output]
data_dir = sys.argv.pop() # path to the dataset
ROOT_DIR = sys.argv.pop() # path to the root folder
uniprotid_list = sys.argv.pop() # path
path_output = os.path.join(ROOT_DIR, rdf_file)
print("%s .... %s ...%s .... %s" % (data_dir,ROOT_DIR,uniprotid_list,path_output))


annotation = milletParser(data_dir,"annotation")
#print(annotation)
annotation2RDF(annotation,path_output)