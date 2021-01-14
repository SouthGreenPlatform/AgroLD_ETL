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
    count = 0
    ttl_handle = open(path_output, "w")
    prot_handle = open(uniprotid_list, 'w')
    prot_buffer = ''
    ttl_handle.write(str(getRDFHeaders()))

    for geneid in annotation:
        annotation_list = annotation[geneid]
        count += 1
        taxon_id = "4543"
        rdf_buffer = ''
        rdf_buffer += "<" + base_resource_uri + "transcript/" + geneid + ">\n"
        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
        rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + geneid + "\" ;\n"
        rdf_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        rdf_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + geneid + "\" ;\n"
        # rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ensembl_transcript_ns + records[0] + ";\n"
        rdf_buffer += "\t" + obo_ns + "RO_0002202" + "\t\t" + res_ns + geneid + " ;\n"
        # rdf_buffer += "\t" + dc_ns + "description" + "\t" + "\"%s" % (
        #     re.sub('\"|\'', '', str(records[2]))) + "\" ;\n"
        #print(len(annotation_list))
        for i in range(len(annotation_list)):
            annotation_dict = annotation[geneid][i]

            if (annotation_dict['type']== 'GO'):
                for go_tuple in annotation_dict['value']:
                    if go_tuple != '':
                        go_list = go_tuple.split(';')
                        if go_list[0] != '':
                            go_id = re.sub(':', '_', go_list[0])
                            go_id = re.sub('\s', '', go_id)
                            if re.match(r"^GO_[0-9]{7}$", go_id):
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + go_id + " ;\n"
                            print(go_id)


                #print(annotation_dict['value'])
            # elif (annotation_dict['type']== 'InterPro'):
            #     print(annotation_dict['value'])
            # elif (annotation_dict['type']== 'KEGG'):
            #     print(annotation_dict['value'])
            # elif (annotation_dict['type']== 'Swiss-Prot'): # InterPro #KEGG
            #     print(annotation_dict['value'])
            # elif (annotation_dict['type'] == 'TrEMBL'):
            #     print(annotation_dict['value'])
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

        rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
        prot_handle.write(prot_buffer)
        ttl_handle.write(rdf_buffer)
    ttl_handle.close()
    prot_handle.close()
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