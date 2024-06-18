import sys
from globalVars import *
from globalVars import base_vocab_ns
from gffParser import *
from utils import *
import pprint
import re
import os


'''
Created on Dec, 2019
The GFF2RDF module is created as part of the AgroLD project.

This module contains Parsers, RDF converters and generic functions for handling genomeHub data
It runs with several files downloaded from MSU annotation project http://rice.plantbiology.msu.edu/

1 - First define your input and output path to filename 
path = ''
path_output = ''
2 - Pass the arguments to genParser function
ds = parseGFF3(path)

3 - run the program 
RDFConverter(ds, path_output)

@author: larmande
'''

__author__ = 'larmande'



def RDFConverter(ds, type_list,output_file):
    '''
    This function converts the GFF3 file to RDF format
    :param ds:
    :param type_list:
    :param output_file:
    :return: rdf file
    '''
    type_list = type_list
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
    exon_list = list()
    taxon_id = "3555"
    source_project = "RefBeet-1.2.2"
    schema_number = "68"
    ssp = ""
    print("************* RDF conversion begins***********\n")
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
    rdf_writer.write(pr + "\t" + "TAIR:" + "<" + 	"http://arabidopsis.org/servlets/TairObject?accession=" + "> .\n")
    rdf_writer.write(pr + "\t" + "KEGG:" + "<" + identifiers_uri + "kegg/> .\n")
    rdf_writer.write(pr + "\t" + "EC:" + "<" + identifiers_uri + "ec-code/> .\n")
    rdf_writer.write(pr + "\t" + "InterPro:" + "<" + identifiers_uri + "interpro/> .\n")
    rdf_writer.write(pr + "\t" + "TrEMBL:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + "MSU:" + "<" + identifiers_uri + "ricegap/> .\n")
    rdf_writer.write(pr + "\t" + "SwissProt:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + faldo_ns + "<" + faldo + "> .\n")
    rdf_writer.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n\n")

    for records in ds:
        line_number += 1
        # loop over the file
        #if records['type'] == "chromosome":
        if not records['seqid'] in chromosome_list:
            genome_buffer = ""
            chromosome_list.append(records['seqid'])
            chromosome_nb = getChromosomeNumber(records['seqid'])
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

        # filtering for gene entries only and avoiding duplicates in the list of genes
        # testing the presence of the gene in the type_list

        if records['type'] == "gene" and "gene" in type_list:
            if not records['attributes']['ID'] in gene_list:
                genome_buffer = ''
                (strand,position) = getStrandValue(records['strand'])
                strand = str(strand)
                # print the corresponding gene associated at mRNAs
                # os_japonica_buffer = ''
                gene_list.append(records['attributes']['ID'])
                genome_buffer += base_resource_ns + records['attributes']['ID'] + "\n"
                genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
                if 'Name' in records['attributes']:
                    genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                if 'Note' in records['attributes']:
                    genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes']['Note'] + "\" ;\n"
                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                                 + chromosome_nb + ":" + \
                                 str(records['start']) + "-" + str(records['end']) + ":" + \
                                 strand + "> ;\n"
                genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
                genome_buffer += getFaldoRegion(taxon_id, ssp, chromosome_nb, records['start'], records['end'],
                                                records['strand'])
                rdf_writer.write(genome_buffer)
        # filtering for mRNA
        if records['type'] == "mRNA" and "mRNA" in type_list:

            if not records['attributes']['ID'] in mRNA_list:
                genome_buffer = ''
                (strand,position) = getStrandValue(records['strand'])
                strand = str(strand)
                go_list = list()
                # print the corresponding gene associated at mRNAs
                mRNA_list.append(records['attributes']['ID'])
                genome_buffer += base_resource_ns + records['attributes']['ID'] + "\n"
                genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
                if 'Name' in records['attributes']:
                    genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                if 'Note' in records['attributes']:
                    genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                        'Note'] + "\" ;\n"

                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                if 'Parent' in records['attributes']:
                    parent = records['attributes']['Parent'].split('.')[0]
                    genome_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t" + base_resource_ns + parent +";\n"
                if 'Dbxref' in records['attributes']:
                    for terms in records['attributes']['Dbxref'].split(','):
                        genome_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" +  terms + " ;\n"
                if 'Ontology_term' in records['attributes']:
                    for terms in records['attributes']['Ontology_term'].split(','):
                        if 'GO:' not in terms:
                            genome_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"" + terms.split(':')[1] + "\" ;\n"
                        else:
                            if terms not in go_list:
                                genome_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + \
                                             re.sub(':', '_', terms) + " ;\n"
                                go_list.append(terms)
                genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                                 + chromosome_nb + ":" + \
                                 str(records['start']) + "-" + str(records['end']) + ":" + \
                                 strand + "> ;\n"
                genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
                genome_buffer += getFaldoRegion(taxon_id, ssp, chromosome_nb,records['start'],records['end'],records['strand'])
                rdf_writer.write(genome_buffer)

        if records['type'] == "polypeptide" and "polypeptide" in type_list:
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
path = '/Users/pierre/Downloads/Beta_vulgaris.RefBeet-1.2.2.57.gff3'
path_output = '/Users/pierre/Downloads/Beta_vulgaris.RefBeet-1.2.2.57.rdf'
ds= parseGFF3(path)
#print(ds)
type_list = ['gene','mRNA','polypeptide']
RDFConverter(ds, type_list,path_output)