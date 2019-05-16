import sys
print sys.path
from globalVars import *
from globalVars import base_vocab_ns
from gffParser import *
import pprint
import re
import os
import pandas as pd
import numpy as np
'''
Created on May, 2017
The msuModel module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling MSU data


@author: larmande
'''
__author__  = "larmande"

#TODO 1) Add documentation
#TODO 2) needs code Refactoring
#TODO 3) better Error handling


# tigr_uri = 'http://www.southgreen.fr/agrold/tigr.locus/'
# tigr_ns = 'tigr:'
#
# tigr_g_uri = 'http://identifiers.org/ricegap/'
# tigr_g_ns = 'tigr_gene:'
def geneParser(infile):
    #    pp = pprint.PrettyPrinter(indent=4)
    gene_hash = {}
    tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.\d$')
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    array = pd.read_csv(infile, sep="\t", delimiter=None, dtype='str')
    #array['locus_id'].replace('', np.nan, inplace=True)
    #array.dropna(subset=['locus_id'], inplace=True)
    print array
    return array

def msuModeleRDF(msu_ds, output_file):
  # The differentes variable declaration
    os_japonica_buffer = ''    # initilised the buffer at zero
    line_number = 0
    rdf_writer = open(output_file, "w")
    gene_list = list()

# The first wrinting in the file is the prefix


    print ("*************msu RDF conversion begins***********\n")
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    rdf_writer.write(pr + "\t" + tigr_g_ns + "<" + tigr_g_uri + "> .\n")
    rdf_writer.write(pr + "\t" + tigr_ns + "<" + tigr_uri + "> .\n")

    # Ajout du prefix pour la realese des donnees
    rdf_writer.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")

# In here we buil the modele and writer in file with ttl format

    for records in msu_ds.as_matrix(columns=None):
        line_number+=1
        # Mapping Triple
        if not (records[0] == "None" or records[1] == "None"):
            print("test passed \n")
            os_japonica_buffer = ''
            os_japonica_buffer += tigr_g_ns + records[1].split('.')[0] + "\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + records[0] + " ;\n"
            # os_japonica_buffer += "\t" + rdfs_ns + "has_dbxref" + "\t" + ensembl_ns + records[0] + " ;\n"
            os_japonica_buffer += "\t" + owl_ns + "sameAs" + "\t" + ensembl_ns + records[0]+ " .\n\n"


            os_japonica_buffer += ensembl_ns + records[0] + "\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + tigr_g_ns + records[1].split('.')[0] + " ;\n"
            # os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + tigr_g_ns + records[1].split('.')[0] + " ;\n"
            os_japonica_buffer += "\t" + owl_ns + "sameAs" + "\t" + tigr_g_ns + records[1].split('.')[0] + " .\n\n"

            #os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
            for gene in records[1].split(','):
                os_japonica_buffer += tigr_ns + gene + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + tigr_g_ns + records[1].split('.')[0] + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + records[0] + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + tigr_g_ns + records[1].split('.')[0] + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + ensembl_ns + records[0] + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + tigr_g_ns + records[1].split('.')[0] + " .\n\n"
            rdf_writer.write(os_japonica_buffer)
            print(os_japonica_buffer)


    print(line_number)


pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
#path = '/Users/plarmande/Downloads/IRGSP-1.0_representative/transcripts_mRNA.gff'
path_output = '/Users/plarmande/Downloads/RAP-MSU_2019-03-22.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
#ds = parseGFF3(path)   # The parsing file withe tropGeneParser()
#pp.pprint(ds)    # For to see in teminal the parsing

#ds = os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)
ds = geneParser('/Users/plarmande/Downloads/RAP-MSU_2019-03-22.txt')
msuModeleRDF(ds, path_output)
