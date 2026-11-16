import sys
import argparse
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
def geneParser(infile, header="infer"):
    """header must be None for regulation_merged_*.txt -- unlike the TF list,
    it has no header row, and its first row is real data that pandas would
    otherwise silently swallow as column names."""
    #    pp = pprint.PrettyPrinter(indent=4)
    tigr = re.compile(tigr_pattern)
    rap = re.compile(rap_pattern)
    array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str', header=header)
    #array['locus_id'].replace('', np.nan, inplace=True)
    #array.dropna(subset=['locus_id'], inplace=True)
    return array

def RDFConverter(ds, output_file, taxon_id):
    os_japonica_buffer = ''  # initilised the buffer at zero
    line_number = 0
    rdf_writer = open(output_file, "w")
    gene_list = list()
    mRNA_list = list()
    fam_list = list()

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
            buffer += "<" + base_resource_uri + "gene_family/" + records[2] + ">\n"
            buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
            buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Transcription_Factor" + " ;\n"
            buffer += "\t" + rdfs_ns + "label" + "\t" +  "\"" + records[2] +"\" ;\n"
            buffer += "\t" + dc_ns + "identifier" + "\t" +  "\"" + records[2] +"\" ;\n"
            # predicate hasMember -- records[0] is the protein/translation id, records[1] the gene id
            buffer += "\t" + obo_ns + "RO_0002351" + "\t" + ensembl_protein_ns + records[0] + " ;\n"
            buffer += "\t" + obo_ns + "RO_0002351" + "\t" + base_resource_ns + records[1] + " ;\n"
            buffer = re.sub(' ;$', ' .\n', buffer)
            fam_list.append(family_name)
        else:
            buffer += "<" + base_resource_uri + "gene_family/" + records[2] + ">" + "\t" + obo_ns + "RO_0002351" + "\t" + \
                      ensembl_protein_ns + records[0] + ".\n"
            buffer += "<" + base_resource_uri + "gene_family/" + records[2] + ">" + "\t" + obo_ns + "RO_0002351" + "\t" + \
                      base_resource_ns + records[1] + " .\n"
                # mRNA uri isMemberOf family uri
        buffer += "<" + ensembl_protein_uri + records[0] + ">"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Transcription_Factor" + " ;\n"
        buffer += "\t" + obo_ns + "RO_0002350" + "\t" + "<" + base_resource_uri + "gene_family/" + records[2] + ">" + " .\n"
        buffer += "<" + base_resource_uri + records[1] + ">"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Transcription_Factor" + " ;\n"
        buffer += "\t" + obo_ns + "RO_0002350" + "\t" + "<" + base_resource_uri + "gene_family/" + records[2] + ">" + " .\n\n"
        # gene uri isMemberOf family uri > inference ?

        buffer = re.sub(' ;$', ' .\n', buffer)

        rdf_writer.write(buffer)
        print(buffer)

    print("*************** PlantTFDB RDF conversion completed ************\n")

def RDFRegulate(ds, regulation_rdf):
    # No taxon_id here: unlike RDFConverter, this function never emitted one
    # (the identically-named local variable it used to declare was dead code).
    buffer = ''  # initilised the buffer at zero
    line_number = 0
    rdf_writer = open(regulation_rdf, "w")
    gene_list = list()
    mRNA_list = list()
    fam_list = list()

    print("************* RDF conversion begins***********\n")
    rdf_writer.write(str(getRDFHeaders()))
    for records in ds.to_numpy():
        buffer = ''
        buffer += "<" + base_resource_uri + records[0] + ">"
        buffer += "\t" + obo_ns + "RO_0002211" + "\t" + "<" + base_resource_uri + records[2] + ">" + " .\n"
        # buffer += "\t" + base_vocab_ns + "regulationType"
#      http://purl.obolibrary.org/obo/RO_0002211
        rdf_writer.write(buffer)
        print(buffer)
def main():
    parser = argparse.ArgumentParser(
        description="Convert PlantTFDB gene family and regulation files to RDF.")
    parser.add_argument("tf_list", help="transcription factor family list (e.g. Mes_TF_list.txt)")
    parser.add_argument("tf_output", help="output Turtle file for the gene family RDF")
    parser.add_argument("regulation_file", help="regulation file (e.g. regulation_merged_Mes.txt)")
    parser.add_argument("regulation_output", help="output Turtle file for the regulation RDF")
    parser.add_argument("--taxon-id", required=True,
                        help="NCBI taxon id of the species this data is for "
                             "(e.g. 39947 for Oryza sativa Japonica) -- the "
                             "input file carries no species information, so "
                             "the caller must supply it")
    args = parser.parse_args()

    ds = geneParser(args.tf_list)
    RDFConverter(ds, args.tf_output, args.taxon_id)
    ds2 = geneParser(args.regulation_file, header=None)
    RDFRegulate(ds2, args.regulation_output)


if __name__ == "__main__":
    main()
