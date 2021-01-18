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

def fileParser(infile,list, mapping, outfile):
    #array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str')
    array1 = pd.read_csv(list, sep="\t", delimiter=None, dtype='str' , header=None, skip_blank_lines=True)
    goarray = pd.read_csv(mapping, sep=" > ", delimiter=None, dtype='str' , header=None, engine='python', skip_blank_lines=True)
    godic = {}
    for line in goarray.to_numpy():
        iprtuple = line[0]
        iprid = iprtuple.split(' ')[0]
        iprid = re.sub('InterPro:', '', iprid)
        item = line[1]
        goitem = item.split('; ')[1]
        if iprid in godic:
            godic[iprid].append(goitem)
        else:
            godic[iprid] = []
            godic[iprid].append(goitem)

    rdf_writer = open(outfile, "w")
    rdf_writer.write(str(getRDFHeaders()))
    nb_line=0
    for records in array1.to_numpy():
        iprid = records[0]
        nb_line +=1
        buffer = ''
        buffer += "<" + interpro_ns + iprid + ">\n"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Interpro_domain" + " ;\n"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + records[1] + " ;\n"
        buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + records[2] + "\" ;\n"
        buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + iprid + "\" ;\n"
        if iprid in godic:
            for goid in godic[iprid]:
                goid = re.sub(':','_',goid)
                if goid != '':
                    buffer += "\t" + base_vocab_ns + "mappingWith" + "\t" + obo_ns + goid + " ;\n"
        buffer = re.sub(' ;$', ' .\n', buffer)
        rdf_writer.write(buffer)
        print(buffer)
    print("nb de lignes traitees:" + str(nb_line))

    with open(infile, 'r', encoding='utf-8') as file:
        root_item = ""
        current_item = ""
        second_item  = ""
        third_item = ""
        fourth_item = ""
        #rdf_writer = open(outfile, "w")
        nb_line = 0
        #rdf_writer.write(str(getRDFHeaders()))
        for line in file:
            buffer = ''
            if not line.startswith('--'):
                nb_line += 1
                key,value = line.split('::')
                root_item = key
            elif not line.startswith('----'):
                nb_line += 1
                key, value = line.split('::')
                key = re.sub('--','',key)
                second_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + root_item + " ;\n"
            elif not line.startswith('------'):
                nb_line += 1
                key, value = line.split('::')
                key = re.sub('--', '', key)
                third_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + second_item + " ;\n"
            elif not line.startswith('--------'):
                nb_line += 1
                key, value = line.split('::')
                key = re.sub('--', '', key)
                fourth_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + third_item + " ;\n"
            else:
                nb_line += 1
                key, value = line.split('::')
                key = re.sub('--', '', key)
                parent_item = current_item
                current_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + fourth_item + " ;\n"

            buffer = re.sub(' ;$', ' .\n', buffer)
            rdf_writer.write(buffer)
            print(buffer)
        print("nb de lignes traitees:" + str(nb_line))



#TEST PARAM
path1 = '/Users/pierre/workspace2015/datasets/ParentChildTreeFile.txt'
path2 ='/Users/pierre/workspace2015/datasets/entry.list'
path3 = '/Users/pierre/workspace2015/datasets/interpro2go'
path_output = '/Users/pierre/workspace2015/datasets/interpro.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
#ds = geneParser(path)   # The parsing file
#pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

#RDFConverter(ds, path_output)

fileParser(path1,path2,path3,path_output)