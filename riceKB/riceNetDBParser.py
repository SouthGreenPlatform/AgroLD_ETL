#!/usr/bin/env python
'''
Created on June 6, 2018
The riceNetDB module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling riceNetDB data

TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: venkatesan
'''
import pprint
from globalVars import *
import re
import os
import sys
#print sys.path
from globalVars import *
from globalVars import base_vocab_ns
import pprint
import re
import os
import pandas as pd
import numpy as np
import urllib
import json
'''
OryzaBase Fields
        gene
        description	
        os	
        chr	
        start	
        stop	
        nl	
        osa	
        dosa	
        NCBI	
        kegg	
        metacyc	
        alternative	
        enzyme	
        enzymedes	
        reaction	
        pubmed	
        uniprot	
        miRNA	
        mirs
'''
def riceNeTDBParser(infile):
    #    pp = pprint.PrettyPrinter(indent=4)
    gene_hash = {}
    tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.\d$')
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    tair_pattern = re.compile(r'^AT[1-5]G\d{5}$')
    prot_pattern = re.compile(
        r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')
    ont_pattern = re.compile(r'^\w+\:\d{7}$')
    array = pd.read_csv(infile, sep=",", dtype='str')
    #array['locus_id'].replace('', np.nan, inplace=True)
    #array.dropna(subset=['locus_id'], inplace=True)
    print array.head()
    return array




def riceNetGeneDBModelRDF(ricenetdb_ds, output_file):
  # The differentes variable declaration
    os_japonica_buffer = ''    # initilised the buffer at zero

    line_number = 0
    rdf_writer = open(output_file, "w")


    pubmed_pattern = re.compile(r'^\d+$')
    ncbi_pattern = re.compile(r'^[A-Z]{2}\d{6}$')

# The first wrinting in the file is the prefix


    print ("*************RiceNetDB RDF conversion begins***********\n")
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    rdf_writer.write(pr + "\t" + mirbase_ns + "<" + mirbase_uri + "> .\n")
    rdf_writer.write(pr + "\t" + mirbase_mature_ns + "<" + mirbase_mature_uri + "> .\n")
    rdf_writer.write(pr + "\t" + protein_ns + "<" + protein_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    rdf_writer.write(pr + "\t" + mRNA_ns + "<" + mRNA_uri + "> .\n")
    rdf_writer.write(pr + "\t" + cDNA_ns + "<" + cDNA_uri + "> .\n")
    rdf_writer.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n")
    rdf_writer.write(pr + "\t" + interpro_ns + "<" + interpro_uri + "> .\n")
    rdf_writer.write(pr + "\t" + embl_ns + "<" + embl_uri + "> .\n")
    rdf_writer.write(pr + "\t" + uniprot_ns + "<" + uniprot_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ncbi_gene_ns + "<" + ncbi_gene_uri + "> .\n")
    rdf_writer.write(pr + "\t" + pubmed_ns + "<" + pubmed_uri + "> .\n")
    rdf_writer.write(pr + "\t" + tenor_ns + "<" + tenor_uri + "> .\n")
    rdf_writer.write(pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n")
    rdf_writer.write(pr + "\t" + tigr_g_ns + "<" + tigr_g_uri  + "> .\n")
    rdf_writer.write(pr + "\t" + RiceNetDB_gene_ns + "<" + RiceNetDB_gene_ns + "> .\n")
    rdf_writer.write(pr + "\t" + RiceNetDB_protein_ns + "<" + RiceNetDB_protein_uri + "> .\n")
    # Ajout du prefix pour la realese des donnees
    rdf_writer.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")

# In here we buil the modele and writer in file with ttl format
#TODO check URI & Namespace LOC_Os01g01010.1 go to http://identifiers.org/ricegap/ or its LOC_Os01g01010
#TODO wich predicate for annotation ? Comment or Description or Annotation
#TODO uniformiser entites choromosomes les valeurs integer ou string ? associer la source.
#TODO has_synonym do we set up URIs for linking ?
    for index, records in ricenetdb_ds.iterrows():
        print records['gene'] + "\n"
        print records['gene'].split(".")[0] + "\n"

        line_number+=1
        os_japonica_buffer = ''
        os_japonica_buffer += RiceNetDB_gene_ns + records['gene'] + "\n"
        os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Gene" + " ;\n"
        os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"

        os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['gene'] + "\" ;\n"
        #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000234" + " ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "description" + "\t" + '"%s"' % (records['description']) + " ;\n"
        if records['os'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_rap_identifier" + "\t" + rapdb_ns + records['os'] + " ;\n"

        os_japonica_buffer += "\t" + base_vocab_ns + "derives_from" + "\t" + tigr_g_ns + records['gene'].split(".")[0] + " ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['stop']) + "\"^^xsd:integer ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "sequence_length" + "\t" + " \"" + str(
            records['nl']) + "\"^^xsd:integer ;\n"

        #os_japonica_buffer += "\t" + base_vocab_ns + "develops_from" + "\t\t" + ensembl_ns + records['attributes']['Locus_id'] + " ;\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + "" + chromosome_ns + re.sub('Chr', '',
                                                                                                            records[
                                                                                                                'chr']) + " ;\n"
        if records['osa'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + records['osa'] + " ;\n"
        if records['dosa'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + records['dosa'] + " ;\n"
        if records['NCBI'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + ncbi_gene_ns + records['NCBI'] + " ;\n"
        if records['protein'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + RiceNetDB_protein_ns + records['protein'] + " ;\n"
        if records['kegg'] is not np.nan:
            for kegg_path_term in records['kegg'].split(";"):
                os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_path_ns+ kegg_path_term + " ;\n" #http://identifiers.org/kegg.pathway/
        if records['metacyc'] is not np.nan:
            for metacyc_term in records['metacyc'].split(";"):
                os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_ns + metacyc_term + " ;\n"
        if records['alternative'] is not np.nan:
            for syn_term in records['alternative'].split(";"):
                #syn_term = re.sub('"', '', syn_term)
                os_japonica_buffer += "\t" + base_vocab_ns + "has_synonym" + "\t" + '"%s"' % (syn_term) + " ;\n"
                #os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_ns + records['alternative'] + " ;\n"
        if records['enzyme'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_ns + records['enzyme'] + " ;\n"
        if records['enzymedes'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_ns + records['enzymedes'] + " ;\n"
        if records['reaction'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + kegg_ns + records['reaction'] + " ;\n"
        if records['pubmed'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + pubmed_ns + records['pubmed'] + " ;\n"
        if records['uniprot'] is not np.nan:
            os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + uniprot_ns + records['uniprot'] + " ;\n"
        if records['miRNA'] is not np.nan:
            for term in records['miRNA'].split(";"):
                os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + mirbase2_ns + term + " ;\n"
        if records['mirs'] is not np.nan:
            for term in records['mirs'].split(";"):
                os_japonica_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + mirbase2_ns + term + " ;\n"

        print(os_japonica_buffer)
        os_japonica_buffer = re.sub(' ;\n$', ' .\n\n', os_japonica_buffer)
        rdf_writer.write(os_japonica_buffer)


    print(line_number)


pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
path = '/Users/plarmande/PycharmProjects/RiceNetDB/gene.csv'
path_output = '/Users/plarmande/PycharmProjects/RiceNetDB/ricenetdb_gene.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
ds = riceNeTDBParser(path)   # The parsing file withe tropGeneParser()
#pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

riceNetGeneDBModelRDF(ds, path_output)
