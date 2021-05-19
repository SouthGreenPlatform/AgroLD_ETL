import sys
print(sys.path)
from riceKB.globalVars import *
from riceKB.globalVars import base_vocab_ns
from riceKB.utils import *
import re
import os
import pandas as pd
import numpy as np


'''
Created on May, 2017
Updated on Dec, 2021
The QtaroParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Qtaro data

TODO:
    1) Add documentation
    2) better Error handling
@author: larmande
'''
__author__  = "larmande"


# ogro_id;gene;gene_symbol;character_major;character_minor;chromosome;start;end;locus_id;browse;isolation;objective;doi

taxon_id = "39947"

def geneParser(infile):
    #    pp = pprint.PrettyPrinter(indent=4)
    gene_hash = {}
    tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.\d$')
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    tair_pattern = re.compile(r'^AT[1-5]G\d{5}$')
    prot_pattern = re.compile(
        r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')
    ont_pattern = re.compile(r'^\w+\:\d{7}$')
    array = pd.read_csv(infile, sep="|", delimiter=None, quotechar = '"', encoding = "ISO-8859-1", dtype='str')
    array['locus_id'].replace('', np.nan, inplace=True)
    array.dropna(subset=['locus_id'], inplace=True)
    print(array)
    return array

#id|qtl_gene_name|character_major|character_minor|marker_volume|chromosome|start|end|source|lod|crossed_a|crossed_b|direction|references|ref_number|marker|marker_physical|marker_fine_1|marker_fine_2|marker_fine_3|marker_interval_1|marker_interval_2|marker_interval_3|marker_single|analytic_group|group_size|trait_name|contribution_rate|synergy_effect|publication_year

def qtlParser(infile):
    raw_headers = "id|qtl_gene_name|character_major|character_minor|marker_volume|chromosome|start|end|source|lod|crossed_a|crossed_b|direction|references|ref_number|marker|marker_physical|marker_fine_1|marker_fine_2|marker_fine_3|marker_interval_1|marker_interval_2|marker_interval_3|marker_single|analytic_group|group_size|trait_name|contribution_rate|synergy_effect|publication_year"
    headers = raw_headers.split("|")
    qtl_ds = list()

    #    pp = pprint.PrettyPrinter(indent=4)
    qtl_ds = pd.read_csv(infile, sep="|", delimiter=None, header=None, dtype='str', encoding = "ISO-8859-1")

    # fileHandle = open(infile, "r")
    # lines = fileHandle.readlines()
    # lines.pop(0)  # remove header
    # for line in lines:
    #     line = re.sub('\n$', '', line)
    #     items = line.split('\t')
    #     qtl_ds.append(dict(zip(headers, items)))
    #
    # fileHandle.close()
    print(qtl_ds)
    return qtl_ds


'''
 RDF Converters
'''
# ogro_id;gene;gene_symbol;character_major;character_minor;chromosome;start;end;locus_id;browse;isolation;objective;doi
def qtaroGeneRDF(infile, output_dir):
    gene_buffer = ''
    to_hash = dict()
    gene_counter = 0
    turtle_file_name = "qtaro.genes.ttl"
    outfile = os.path.join(output_dir, turtle_file_name)
    outHandle = open(outfile, "w")
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    taxon_id = "39947"
    print("*********** Parsing Qtaro Gene data ***************\n")

    gene_ds = geneParser(infile)

    print("************* Qtaro Gene RDF conversion begins***********\n")

    '''
    Ajout du prefix pour la release des donnees
    '''
    # outHandle.write(base + "\t" + "<" + base_uri + "> .\n")
    # outHandle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    # outHandle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    # outHandle.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    # outHandle.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    # outHandle.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    # outHandle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    # outHandle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    # outHandle.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    # outHandle.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    # outHandle.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n\n")
    outHandle.write(str(getRDFHeaders()))
    for records in gene_ds.to_numpy():
        gene_buffer = ''
        gene_counter += 1

        # ogro_id;gene;gene_symbol;character_major;character_minor;chromosome;start;end;locus_id;browse;isolation;objective;doi
        print(records)
        if rap_pattern.match(records[8]) and isinstance(records[8],str):
            gene_buffer += base_resource_ns + records[8] + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            label = re.sub('\"+', '', records[1])
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (label) + " ;\n"
            gene_buffer += "\t" + skos_ns + "prefSymbol" + "\t" + '"%s"' % (records[2]) + " ;\n"
            # gene_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t" + '"%s"' % (records[5]) + " ;\n"
            # gene_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + '"%s"' % (records[6]) + " ;\n"
            # gene_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + '"%s"' % (records[7]) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % (records[3]) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % (records[4]) + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + records[8] + " ;\n"
            description = re.sub('\"+', '', records[11])
            gene_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % (description) + " ;\n"
            if records[12] is not np.nan:
                gene_buffer += "\t" + dc_ns + "references" + "\t" +  "<" + doi_uri + records[12] + "> ;\n"
            gene_buffer = re.sub(' ;$', ' .\n', gene_buffer)
            outHandle.write(gene_buffer)
    outHandle.close()
    print("Number of Genes: %s\n" % (str(gene_counter)))
    print("********* Qtaro GENE RDF completed ***********\n")

def qtaroQTLRDF(infile, output_dir):
    qtl_buffer = ''
    to_hash = dict()
    qtl_counter = 0
    turtle_file_name = "qtaro.qtl.ttl"
    outfile = os.path.join(output_dir, turtle_file_name)
    outHandle = open(outfile, "w")

    print("*********** Parsing Qtaro QTL data ***************\n")

    qtl_ds = qtlParser(infile)


    #    print "Gramene QTL data has been parsed!\n"
    #    print "*************************************\n"

    print("************* Qtaro QTL RDF conversion begins***********\n")

    outHandle.write(base + "\t" + "<" + base_uri + "> .\n")
    outHandle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    outHandle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    outHandle.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    outHandle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    outHandle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    outHandle.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")

    '''
    Ajout du prefix pour la release des donnees
    '''
    outHandle.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n\n")

    for records in qtl_ds.to_numpy():
        qtl_buffer = ''
        qtl_counter += 1
        #chrm = records['Chromosome'].replace("Chr. ", "")
        #to_id = records['TOid'].replace(":", "_")
        trait1 = re.sub('"', '',records[2])
        trait2 = re.sub('"', '',records[3])
        trait3 = re.sub('"', '', str(records[26]))
        print(records)
        qtl_buffer += "<" + base_resource_uri + "qtaro.qtl/" + records[0] + ">\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "QTL" + " ;\n"
        # qtl_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
        # qtl_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + qtl_term + " ;\n"
        qtl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (records[1]) + " ;\n"
        # URI du Chromosome ensembl http://www.southgreen.fr/agrold/resource/oryza_sativa/IRGSP-1.0/4
        #
        qtl_buffer += "\t" + faldo_ns + "location" + "\t" + "<"+ chromosome_uri + taxon_id + "/{}:{}-{}:1>".format(records[5],records[6],records[7]) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % (trait1) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % (trait2) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % (trait3) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "lod" + "\t" + '"%s"' % (records[9]) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "date" + "\t" + '"%s"' % (records[29]) + " .\n\n"
        # Region
        qtl_buffer += "<"+ chromosome_uri + taxon_id + "/{}:{}-{}:1>".format(records[5],records[6],records[7])  + "\n"
        qtl_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + taxon_id + "/{}:{}-{}:1>".format(records[5],records[6],records[7])  + "\" ;\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
        qtl_buffer += "\t" + faldo_ns + "begin" + "\t" +  "<"+ chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5],records[6])  + " ;\n"
        qtl_buffer += "\t" + faldo_ns + "end" + "\t" + "<"+ chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5],records[7]) + "  .\n\n"

        # Position 1
        qtl_buffer += "<"+ chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5],records[6])  + "\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ForwardStrandPosition" + " ;\n"
        qtl_buffer += "\t" + faldo_ns + "position" + "\t" + records[6] + " ;\n"
        qtl_buffer +=   "\t" + faldo_ns + "reference" + "\t" + "<"+ chromosome_uri + taxon_id + "/{}>".format(records[5])  + ". \n\n"
        # Position 2
        qtl_buffer += "<" + chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5], records[7]) + "\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ForwardStrandPosition" + " ;\n"
        qtl_buffer += "\t" + faldo_ns + "position" + "\t" + records[7] + " ;\n"
        qtl_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + "/{}>".format(
            records[5]) + ". \n\n"

        #if records['marker_interval_1']

        #        if to_id not in to_hash:
        #            outHandle.write(obo_ns + to_id + "\n")
        #            outHandle.write("\t" + rdf_ns + "type" + "\t" + obo_ns + plant_trait_term + " ;\n") #base_vocab_ns + "Concept"
        #            outHandle.write("\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + plant_trait_term + " ;\n")
        #            outHandle.write("\t" + rdfs_ns + "label" + "\t" + '"%s"' % (records['TraitName']) + " ;\n")
        #            outHandle.write("\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['TraitSymbol']) + " ;\n")
        #            outHandle.write("\t" + base_vocab_ns + "has_category" + "\t" + '"%s"' % (records['Category']) + " .\n")
        #            to_hash[to_id] = 1
        qtl_buffer = re.sub(' ;$', ' .\n', qtl_buffer)
        outHandle.write(qtl_buffer)
    outHandle.close()
    print("Number of QTLs: %s\n" % (str(qtl_counter)))
    print("********* Qtaro QTL RDF completed ***********\n")

#geneParser('../test_files/qtaro/Qtaro-Gene-export.csv')

#qtaroGeneRDF('/Users/pierre/workspace2015/datasets/qtaro_gene.csv','/Users/pierre/workspace2015/datasets')
qtaroQTLRDF('/Users/pierre/workspace2015/datasets/qtaro.csv','/Users/pierre/workspace2015/datasets')