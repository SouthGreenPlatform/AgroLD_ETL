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
Created on June, 2021
This module is created as part of the AgroLD Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Ensembl genome annotation \
datasets to normalize with AgroLD

It runs with several files downloaded from Ensembl Plants project available at  ftp.ensemblgenomes.org


@author: larmande
'''


__author__  = "larmande"
def printHeader(rdf_writer):
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl_uri + "> .\n")
    rdf_writer.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    rdf_writer.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    rdf_writer.write(pr + "\t" + protein_ns + "<" + protein_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    rdf_writer.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n")
    rdf_writer.write(pr + "\t" + interpro_ns + "<" + interpro_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n")
    rdf_writer.write(pr + "\t" + "Phytozome:" + "<" + identifiers_uri + "phytosome/> .\n")
    rdf_writer.write(pr + "\t" + "PFAM:" + "<" + identifiers_uri + "pfam/> .\n")
    rdf_writer.write(pr + "\t" + "Panther:" + "<" + identifiers_uri + "panther.family/> .\n")
    rdf_writer.write(pr + "\t" + "KOG:" + "<" + identifiers_uri + "kog/> .\n")
    rdf_writer.write(pr + "\t" + "TAIR:" + "<" + "http://arabidopsis.org/servlets/TairObject?accession=" + "> .\n")
    rdf_writer.write(pr + "\t" + "KEGG:" + "<" + identifiers_uri + "kegg/> .\n")
    rdf_writer.write(pr + "\t" + "EC:" + "<" + identifiers_uri + "ec-code/> .\n")
    rdf_writer.write(pr + "\t" + "InterPro:" + "<" + identifiers_uri + "interpro/> .\n")
    rdf_writer.write(pr + "\t" + "TrEMBL:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + "MSU:" + "<" + identifiers_uri + "ricegap/> .\n")
    rdf_writer.write(pr + "\t" + "SwissProt:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + faldo_ns + "<" + faldo + "> .\n")
    # Ajout du prefix pour la release des donnees
    rdf_writer.write(pr + "\t" + "gene:" + "<" + base_resource_uri + "> .\n")
    rdf_writer.write(pr + "\t" + "transcript:" + "<" + base_resource_uri + "transcript/" + "> .\n")
    rdf_writer.write(pr + "\t" + "CDS:" + "<" + base_resource_uri + "CDS/" + "> .\n")
    rdf_writer.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n")

def geneEntityWriter(type,records,gene_list,taxon_id,chromosome_nb,ssp, source_project):
    if not records['attributes']['gene_id'] in gene_list:
        genome_buffer = ''
        (strand, position) = getStrandValue(records['strand'])
        strand = str(strand)
        # print the corresponding gene associated at mRNAs
        # os_japonica_buffer = ''
        gene_list.append(records['attributes']['gene_id'])
        # genome_buffer += base_resource_ns + records['attributes']['gene_id'] + "\n"
        genome_buffer += records['attributes']['ID'] + "\n"
        genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
        genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
        if 'Name' in records['attributes']:
            genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
        else:
            genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['gene_id'] + "\" ;\n"
        if 'description' in records['attributes']:
            genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                'description'] + "\" ;\n"
        genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id + "/" + ssp \
                         + chromosome_nb + ":" + \
                         str(records['start']) + "-" + str(records['end']) + ":" + \
                         strand + "> ;\n"
        genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
        genome_buffer += getFaldoRegion(taxon_id, chromosome_nb, records['start'], records['end'],
                                        records['strand'])
        return genome_buffer
def RNAEntityWriter(type,records,mRNA_list,taxon_id,chromosome_nb,ssp, source_project):
    if not records['attributes']['transcript_id'] in mRNA_list:
        genome_buffer = ''
        (strand, position) = getStrandValue(records['strand'])
        strand = str(strand)
        go_list = list()
        # print the corresponding gene associated at mRNAs
        mRNA_list.append(records['attributes']['transcript_id'])
        # genome_buffer += base_resource_ns + "transcript/"+ records['attributes']['transcript_id'] + "\n"
        genome_buffer += records['attributes']['ID'] + "\n"
        genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
        genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
        if 'Name' in records['attributes']:
            genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
        else:
            genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes'][
                'transcript_id'] + "\" ;\n"
        if 'description' in records['attributes']:
            genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                'description'] + "\" ;\n"

        genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        if 'Parent' in records['attributes']:
            parent = records['attributes']['Parent'].split('.')[0]
            genome_buffer += "\t" + sio_ns + 'SIO_010081' + "\t" + parent + ";\n"
            genome_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t" + parent + ";\n"
        if 'Dbxref' in records['attributes']:
            for terms in records['attributes']['Dbxref'].split(','):
                genome_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + terms + " ;\n"
        if 'Ontology_term' in records['attributes']:
            for terms in records['attributes']['Ontology_term'].split(','):
                if 'GO:' not in terms:
                    genome_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"" + terms.split(':')[
                        1] + "\" ;\n"
                else:
                    if terms not in go_list:
                        genome_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + \
                                         re.sub(':', '_', terms) + " ;\n"
                        go_list.append(terms)
        genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id + "/" + ssp \
                         + chromosome_nb + ":" + \
                         str(records['start']) + "-" + str(records['end']) + ":" + \
                         strand + "> ;\n"
        genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
        genome_buffer += getFaldoRegion(taxon_id, chromosome_nb, records['start'], records['end'], records['strand'])
        return genome_buffer
def RDFConverter(ds, output_file):
    os_japonica_buffer = ''  # initilised the buffer at zero
    number_match_part_sbgi = 0
    number_match_part_kome = 0
    number_exon = 0
    number_cds = 0
    line_number = 0
    number_five_prime_UTR = 0
    number_three_prime_UTR = 0
    rdf_writer = open(output_file, "w")
    chromosome_list = list()
    gene_list = list()
    mRNA_list = list()
    taxon_id = "39947"
    source_project = ""
    schema_number = ""
    ssp = ""
    print("************* RDF conversion begins***********\n")
    printHeader(rdf_writer)

    for records in ds:
        line_number += 1
        # loop over the file
        #if records['type'] == "chromosome":
        if not records['seqid'] in chromosome_list:
            genome_buffer = ""
            chromosome_list.append(records['seqid'])
            chromosome_nb = getChromosomeNumber(records['seqid'])
            source_project = records['source']
            if 'Chr' in records['seqid']:
                genome_buffer += "<" + chromosome_uri + taxon_id +"/"+ ssp + chromosome_nb+ ">\n"
                genome_buffer += "\t" +  obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inAssembly" + "\t" + "\"Reference-" + source_project + "\" ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inSchemaNumber" + "\t" + "\""+ schema_number + "\"" + " ;\n"
                genome_buffer = re.sub(' ;$', ' .\n\n', genome_buffer)
                rdf_writer.write(genome_buffer)
            else:
                genome_buffer += "<" + chromosome_uri + taxon_id + "/" + ssp + chromosome_nb + ">\n"
                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inAssembly" + "\t" + "\"Reference-" + source_project + "\" ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inSchemaNumber" + "\t" + "\"" + schema_number + "\"" + " ;\n"
                genome_buffer = re.sub(' ;$', ' .\n\n', genome_buffer)
                rdf_writer.write(genome_buffer)
            print(genome_buffer)
        # filtering for gene entries
        if records['type'] == "gene":
            genome_buffer = geneEntityWriter(records['type'], records, gene_list, taxon_id, chromosome_nb, ssp, source_project)
            rdf_writer.write(genome_buffer)
                # print(genome_buffer)
        # filtering for mRNA
        if records['type'] == "mRNA":
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                             source_project)
            rdf_writer.write(genome_buffer)
                # print(genome_buffer)
        # if records['type'] == "three_prime_UTR":
        #     print('three_prime_UTR')
        # if records['type'] == "five_prime_UTR":
        #     print('five_prime_UTR')
        # if records['type'] == "CDS":
        #     print('CDS')
        # if records['type'] == "exon":
        #     print('exon')
        if records['type'] == "lnc_RNA":
            # print('lnc_RNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "ncRNA_gene":
            # print('ncRNA_gene')
            genome_buffer = geneEntityWriter(records['type'], records, gene_list, taxon_id, chromosome_nb, ssp,
                                             source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "pre_miRNA":
            # print('pre_miRNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "rRNA":
            # print('rRNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "snRNA":
            # print('snRNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "snoRNA":
            # print('snoRNA')
            # print('snRNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        if records['type'] == "tRNA":
            # print('tRNA')
            genome_buffer = RNAEntityWriter(records['type'], records, mRNA_list, taxon_id, chromosome_nb, ssp,
                                            source_project)
            rdf_writer.write(genome_buffer)
        # if records['type'] == "RNase_MRP_RNA":
        #     # print('RNase_MRP_RNA')
        # if records['type'] == "SRP_RNA":
        #     # print('SRP_RNA')
        if records['type'] == "polypeptide":
            genome_buffer = ''
            # print the corresponding gene associated at mRNAs
            # os_japonica_buffer = ''
            gene_list.append(records['attributes']['ID'])
            genome_buffer += base_resource_ns + records['attributes']['ID'] + "\n"
            genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
            genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n"
            if 'Name' in records['attributes']:
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
            else:
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['ID'] + "\" ;\n"
            if 'Note' in records['attributes']:
                genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                    'Note'] + "\" ;\n"
            if 'Derives_from' in records['attributes']:
                genome_buffer += "\t" + base_vocab_ns + "derivesFrom" + "\t"  + base_resource_ns + records['attributes']['Derives_from'] + " ;\n"
            genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
            genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
            rdf_writer.write(genome_buffer)
        # break

pp = pprint.PrettyPrinter(indent=4)




#TEST PARAM
path = '/Users/pierre/workspace2015/datasets/Oryza_sativa.IRGSP-1.0.51.chromosome.1.gff3'
path_output = '/Users/pierre/workspace2015/datasets/Oryza_sativa.IRGSP-1.0.51.chromosome.1.gff3.ttl'
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
ds = parseGFF3(path)   # The parsing file
print(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

RDFConverter(ds, path_output)
