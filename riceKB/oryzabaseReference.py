import gffParser
import pprint
from globalVars import *
from globalVars import base_vocab_ns
import re
import requests
import os, sys
import datetime
import json
import copy
import time
import pandas as pd
import numpy as np
from pandas.io.json import json_normalize
import rdflib
from rdflib.graph import Graph
from rdflib import URIRef
import urllib3
import requests


#pp = pprint.PrettyPrinter(indent=4)

def pubParser(infile):
    pub_hash = {}
    file_reader = open(infile, "r")
    lines = file_reader.readlines()
    lines.pop(0) # remove header
    for line in lines:
        line = re.sub('\n$', '', line)
        records = line.split('\t')
        reference_id = records.pop(0)
        pub_hash[reference_id]= {
                                  'PubMedId': records[0],
                                  'Author': records[1],
                                  'Title': records[2],
                                  'Journal': records[3],
                                  'Volume': records[4],
                                  'Pages': records[5],
                                  'Year': records[6],
                                  'Symbol': [],
                                  'Synonym': []
                                  }
        if records[7]:
                pub_hash[reference_id]['Symbol'] = records[7]
        if records[8]:
            if prot_pattern.match(records[8]):
                pub_hash[reference_id]['ProtID'] = records[8]
    return pub_hash

def printGenes(file):
    pubParser(file)

def pubParserPandas(infile):
    array = pd.read_csv(infile, sep="\t", delimiter=None, dtype='str')
    print array
    return array

def referenceRDF(file, rdf_file, output_dir,type='run'):
    rdf_buffer = ''
    turtle_file = "reference_oryzabase.ttl"
    output_file = os.path.join(output_dir, turtle_file)
    print "*************Parsing %s genome data ***********\n" % (file)

    output_opener = open(output_file, "w+")
    # Printing Prefixes
    output_opener.write(base + "\t" + "<" + base_uri + "> .\n")
    output_opener.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    output_opener.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    output_opener.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    output_opener.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    output_opener.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    output_opener.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    output_opener.write(pr + "\t" + rapdb_gene_ns + "<" + rapdb_gene_uri + "> .\n")
    output_opener.write(pr + "\t" + msu_ns + "<" + msu_uri + "> .\n")
    output_opener.write(pr + "\t" + tair_l_ns + "<" + tair_l_uri + "> .\n")
    output_opener.write(pr + "\t" + up_ns + "<" + uniprot + "> .\n")
    output_opener.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n")
    output_opener.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    output_opener.write(pr + "\t" + faldo_ns + "<" + faldo + "> .\n")
    output_opener.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    output_opener.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    output_opener.write(pr + "\t" + sio_ns + "<" + sio_uri + "> .\n")
    output_opener.write(pr + "\t" + uniprot_ns + "<" + uniprot_uri + "> .\n")
    output_opener.write(pr + "\t" + pubmed_ns + "<" + pubmed_uri + "> .\n\n")

    pub_ds = pubParser(file)
    # global RDFgraph
    # if rdf_file:
    #     RDFgraph = Graph()
    #     RDFgraph.parse(rdf_file, format='turtle')
    #     print('RDF file exists ****')

    print(str(len(pub_ds.keys())) + " Reference parsed")
    print "************* %s RDF conversion begins***********\n" % (file)
    for pub_id in pub_ds:
        if pubmed_pattern.match(pub_ds[pub_id]['PubMedId']):
            print(pub_ds[pub_id]['PubMedId'])
            rdf_buffer += pubmed_ns + pub_id + "\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"


file = '/Users/plarmande/Downloads/Reference_20190617000959.txt'
rdf_file = '/Users/plarmande/Downloads/OryzabaseGeneListEn_20190528010057.ttl'
output_dir = '/Users/plarmande/Downloads/'

referenceRDF(file,rdf_file, output_dir)
#pubParserPandas(file)

#path = '/media/elhassouni/donnees/Noeud-plante-projet/workspace/AgroLD/AgroLD_ETL/test_files/urgi/pseudomolecul_wheat.gff'

#pp.pprint(parseGFF3(path))


