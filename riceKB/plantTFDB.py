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
The genomeHub module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling genomeHub data
It runs with several files downloaded from MSU annotation project http://rice.plantbiology.msu.edu/



3 - run the program 
msuModeleRDF(ds, path_output)

@author: larmande
'''
__author__  = "larmande"

# TODO better Error handling
# TODO modify the help
def geneParser(infile):
    #    pp = pprint.PrettyPrinter(indent=4)
    tigr = re.compile(tigr_pattern)
    rap = re.compile(rap_pattern)
    array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str')
    #array['locus_id'].replace('', np.nan, inplace=True)
    #array.dropna(subset=['locus_id'], inplace=True)
    return array

def RDFConverter(ds, output_file):
    os_japonica_buffer = ''  # initilised the buffer at zero
    line_number = 0
    rdf_writer = open(output_file, "w")
    gene_list = list()
    mRNA_list = list()
    fam_list = list()
    taxon_id = "39947"

    print("************* RDF conversion begins***********\n")
    rdf_writer.write(str(getRDFHeaders()))
    for records in ds.to_numpy():
        line_number += 1
        #print(str(records[2]))
        #for fam in records[2]:
        family_name = re.sub('\n','',records[2]) # re.sub('"', '', records['attributes']['CGSNL Gene Name'])
        buffer = ''
        if family_name not in fam_list:
            # buffer = ''
            buffer += "<" + base_resource_uri + "family/" + records[2] + ">\n"
            buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
            buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Transcription_Factor" + " ;\n"
            buffer += "\t" + rdfs_ns + "label" + "\t" +  "\"" + records[2] +"\" ;\n"
            buffer += "\t" + dc_ns + "identifier" + "\t" +  "\"" + records[2] +"\" ;\n"
            # predicate hasMember
            buffer += "\t" + obo_ns + "RO_0002351" + "\t" + base_resource_ns + records[0] + " ;\n"
            buffer += "\t" + obo_ns + "RO_0002351" + "\t" + base_resource_ns + records[1] + " ;\n"
            buffer = re.sub(' ;$', ' .\n', buffer)
            fam_list.append(family_name)
        else:
            buffer += "<" + base_resource_uri + "family/" + records[2] + ">" + "\t" + obo_ns + "RO_0002351" + "\t" + \
                      base_resource_ns + records[0] + ".\n"
            buffer += "<" + base_resource_uri + "family/" + records[2] + ">" + "\t" + obo_ns + "RO_0002351" + "\t" + \
                      base_resource_ns + records[1] + " .\n"
                # mRNA uri isMemberOf family uri
        buffer += "<" + base_resource_uri + records[0] + ">"
        buffer += "\t" + obo_ns + "RO_0002350" + "\t" + "<" + base_resource_uri + "family/" + records[2] + ">" + " .\n"
        buffer += "<" + base_resource_uri + records[1] + ">"
        buffer += "\t" + obo_ns + "RO_0002350" + "\t" + "<" + base_resource_uri + "family/" + records[2] + ">" + " .\n\n"
        # gene uri isMemberOf family uri > inference ?

        buffer = re.sub(' ;$', ' .\n', buffer)

        rdf_writer.write(buffer)
        print(buffer)

    print("*************** PlantTFDB RDF conversion completed ************\n")


pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
path = '/Users/pierre/workspace2015/datasets/Osj_TF_list.txt'
path_output = '/Users/pierre/workspace2015/datasets/Osj_TF_list.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
ds = geneParser(path)   # The parsing file
pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

RDFConverter(ds, path_output)
