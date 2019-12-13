import sys
print(sys.path)
from riceKB.globalVars import *
from riceKB.globalVars import base_vocab_ns
from riceKB.gffParser import *
import pprint
import re
import os

'''
Created on Dec, 2019
The genomeHub module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling genomeHub data
It runs with several files downloaded from MSU annotation project http://rice.plantbiology.msu.edu/

1 - First define your output path to filename 
path_output = 'all.locus_brief_info.7.ttl'
2 - Pass the arguments to genParser function
ds = geneParser('all.locus_brief_info.7.0.txt',\
                'all.interpro.txt',\
                'all.pfam.txt',\
                'all.GOSlim_assignment.txt')

3 - run the program 
msuModeleRDF(ds, path_output)

@author: larmande
'''
__author__  = "larmande"

#TODO better Error handling

def getStrandValue(strandVar):
    if strandVar == "-":
        strandVar = "-1"
        positionVar = "ReverseStrandPosition"
    else:
        strandVar = "1"
        positionVar = "ForwardStrandPosition"
    return (strandVar,positionVar)

def getFaldoRegion(taxon_id,seqid,start,end,strand):
    (strand, position) = getStrandValue(strand)
    strand = str(strand)
    genome_buffer =''
    # Region

    genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(start) + "-" + str(end) + ":" \
                     + strand + "> \n"
    genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(start) + "-" + str(end) + ":" \
                     + strand + "\";\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
    genome_buffer += "\t" + faldo_ns + "begin" + "\t" + "<" + chromosome_uri + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(start) + ":" + strand + ">;\n"
    genome_buffer += "\t" + faldo_ns + "end" + "\t" + "<" + chromosome_uri + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(end) + ":" + strand + ">  .\n\n"

    # Position 1
    genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(start) + ":" + strand + ">\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
    genome_buffer += "  ;\n"
    genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(start) + " ;\n"
    genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + \
                     "/" + re.sub('Os|Chr', '', seqid) + "> .\n\n"

    # Position 2
    genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
                     + re.sub('Os|Chr', '', seqid) + ":" + \
                     str(start) + ":" + strand + "> \n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
    genome_buffer += " ;\n"
    genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(end) + " ;\n"
    genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + \
                     "/" + re.sub('Os|Chr', '', seqid) + "> .\n\n"
    return genome_buffer