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

'''
Created on May, 2017
The rapdbParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling rapdb data

TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
    4) put @en at label and all string literal
@author: larmande
'''
#  check  agrold_vocabulary:go_term property mapping rdfs:seeAlso + other predicate = en fait uniprot a son propre predicat
#  put @en at label and all string literal - fait pour les labels mais pas pour les autres

__author__  = "larmande"

cds_hash = dict()


def getStrandValue(strandVar):
    if strandVar == "-":
        strandVar = "-1"
        positionVar = "ReverseStrandPosition"
    else:
        strandVar = "1"
        positionVar = "ForwardStrandPosition"
    return (strandVar,positionVar)

def getCDSNumber(transcript, start, number_cds):
    #transcript = re.sub('-', '', transcript)
    if transcript not in cds_hash:
     cds_hash[transcript] = [start]
     number_cds = 1
     return number_cds
    elif start not in cds_hash[transcript]:
        cds_hash[transcript].append(start)
        number_cds += 1
        return number_cds


def rapdbModeleRDF(rapdb_ds, output_file):
  # The differentes variable declaration
    os_japonica_buffer = ''    # initilised the buffer at zero
    number_match_part_sbgi = 0
    number_match_part_kome = 0
    number_exon = 0
    number_cds = 0
    line_number = 0
    ch_number = 0
    chromosome_dict = {}
    number_five_prime_UTR = 0
    number_three_prime_UTR = 0
    rdf_writer = open(output_file, "w")
    chromosome_list = {}
    taxon_id = "39947"



    pubmed_pattern = re.compile(r'^\d+$')
    ncbi_pattern = re.compile(r'^[A-Z]{2}\d{6}$')

    chromosome_size = [43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856]
    # chromosome_size.reverse()
    # The first wrinting in the file is the prefix


    print ("*************RapDB RDF conversion begins***********\n")


    rdf_writer.write(str(getRDFHeaders()))
# In here we buil the modele and writer in file with ttl format

  # Species triple


    os_japonica_buffer = ''
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "OBI_0100026" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + skos_ns + "prefLabel" + "\t\t" + "\"Oryza sativa Japonica Group\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + rdfs_ns + "label" + "\t\t" + "\"Oryza sativa Japonica Group\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + skos_ns + "altLabel" + "\t\t" + "\"Japanese rice\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + dcterms_ns + "identifier" + "\t\t" + taxon_id + " .\n\n"
    # os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + base_vocab_ns + "taxon" + "\t\t" + ncbi_tax_ns + "39947" + " .\n\n"
    print(os_japonica_buffer)
    rdf_writer.write(os_japonica_buffer)
# Chromosome triple
# chromosome:IRGSP-1.0:1:1:43270923:1
#  rdfs:subClassOf obo:SO_0000340
# rdfs:label "Oryza sativa Japonica Group chromosome:IRGSP-1.0:1:1:43270923:1 (IRGSP-1.0)" .
# dc:identifier "chromosome:IRGSP-1.0:1:1:43270923:1" .
# term:inEnsemblAssembly "IRGSP-1.0" .
    for records in rapdb_ds:
        line_number+=1
        if not records['seqid'] in chromosome_list.keys():
            # for chromosome in chromosome_size:
            # inverselist = chromosome_size.reverse()
            chromosome = chromosome_size.pop(0)
            os_japonica_buffer = ''
            ch_number += 1
            print(records['seqid'])
            chromosome_dict[records['seqid']] = {}
            chromosome_dict[records['seqid']] = { 'uri': "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1", 'seqid': records['seqid'], 'number' : str(ch_number), 'nucleotide' : str(chromosome), 'assembly':'IRGSP-1.0'}
            chromosome_list[records['seqid']]= "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1"
            os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1" + "\n"
            os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000340" + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + "Oryza sativa Japonica Group chromosome:" + "IRGSP-1.0:" + str(
                ch_number) + ":1-" + str(chromosome) + ":1" + " (IRGSP-1.0)" + "\"@en ;\n"
            os_japonica_buffer += "\t" + dcterms_ns + "identifier " + "\t" + " \"" + "IRGSP-1.0:" + str(
                ch_number) + ":1-" + str(chromosome) + ":1" + "\" ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "chromosomeNumber" + "\t" + "\"" + str(ch_number) + "\"^^xsd:integer ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "genomeAssembly " + "\t" + " \"" + "IRGSP-1.0" + "\" .\n\n"

            print(os_japonica_buffer)
            # os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
            rdf_writer.write(os_japonica_buffer)

        #RAPDB.gff3

        if records['source']:
            if records['type'] == "gene":
                (strand,position) = getStrandValue(records['strand'])

                os_japonica_buffer = ''

                os_japonica_buffer += rapdb_gene_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" +  ensembl_gene_ns + records['attributes']['ID'] + " ;\n"
                os_japonica_buffer += "\t" + owl_ns + "sameAs" + "\t" + ensembl_gene_ns + records['attributes']['ID'] + " ;\n"

                # rapdb..ID  skos:closeMatch ensembl:id  ## important ##
                os_japonica_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + records['source'] + "\" ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000704" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "hasBiotype" + "\t" + "\"protein_coding\" ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\"@en ;\n"
                os_japonica_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                if 'Note' in records['attributes']:
                    os_japonica_buffer += "\t" + dc_ns + "description" + "\t" + " \"" + records['attributes']['Note'] + "\" ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"

                os_japonica_buffer +=  "\t" + faldo_ns + "location" + "\t"  + chromosome_ns + "IRGSP-1.0:"+ \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:"+ \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(records['start']) + '-' + str(records['end']) + ":" + strand +  "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t"  + " \"" + chromosome_ns + "IRGSP-1.0:"+ \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" +  "\t" + chromosome_ns + "IRGSP-1.0:"+ \
                                      chromosome_dict[records['seqid']]['number']+":"+str(records['start'])+":"+ strand + "  ;\n"
                os_japonica_buffer +=  "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(records['end']) + ":" + strand + "  .\n\n"


                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer +=  "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" +"\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)


            if records['type'] == "mRNA":

                (strand,position) = getStrandValue(records['strand'])
                os_japonica_buffer = ''
                os_japonica_buffer += rapdb_mrna_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000234" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\"@en ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + records['source'] + "\" ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "hasBiotype" + "\t" + "\"protein_coding\" ;\n"
                os_japonica_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ensembl_transcript_ns + records['attributes']['Name'] + ";\n"

                if 'Note' in records['attributes']:
                    os_japonica_buffer += "\t" +  dc_ns + "description" + "\t" + "\"%s" % (records['attributes']['Note']) + "\" ;\n"

                os_japonica_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"
                os_japonica_buffer += "\t" + sio_ns +  "SIO_010081" +  "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"
                # <http://rdf.ebi.ac.uk/resource/ensembl.transcript/Os09t0372700-01> obo:SO_translates_to <http://rdf.ebi.ac.uk/resource/ensembl.protein/Os09t0372700-01> .
                ## os_japonica_buffer += "\t" + obo_ns + "SO_translates_to" + "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] +\
                                      ':' + str(records['start']) + '-' + str(records['end']) + ":" + strand + " ;\n"

                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + "" + chromosome_ns + re.sub('Os', '', records['seqid']) + " ;\n"


                if 'GO' in records['attributes']:
                        for go_term in re.findall(r'GO:[0-9]{7}',records['attributes']['GO']):
                            os_japonica_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + re.sub(':', '_', go_term) + " ;\n"
                            #os_japonica_buffer += "\t" + base_vocab_ns + "comment" + "\t" + '"%s"' % (records['attributes']['GO']) + " ;\n"
                if 'InterPro' in records['attributes']:
                    for ipr_term in re.findall(r'IPR[0-9]{6}', records['attributes']['InterPro']):
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" +"\t" + interpro_ns + ipr_term + " ;\n"
                if 'CGSNL Gene Name' in records['attributes']:
                    syn_term = re.sub('"', '', records['attributes']['CGSNL Gene Name'])
                    os_japonica_buffer += "\t" + skos_ns + "prefLabel" + "\t" + "\"%s" % (syn_term) + "\" ;\n"
                if 'CGSNL Gene Symbol' in records['attributes']:
                    sym_term = re.sub('"', '', records['attributes']['CGSNL Gene Symbol'])
                    os_japonica_buffer += "\t" + skos_ns + "prefSymbol" + "\t" +  "\"%s" % (sym_term) + "\" ;\n"
                if 'Literature_PMID' in records['attributes']:
                    if pubmed_pattern.match(records['attributes']['Literature_PMID']):
                        os_japonica_buffer += "\t" + dc_ns + "references" + "\t\t" + pubmed_ns + records['attributes']['Literature_PMID'] + " ;\n"
                if 'ORF_evidence' in records['attributes']:
                    if '(UniProt)' in records['attributes']['ORF_evidence']:
                        uni_term = records['attributes']['ORF_evidence'].split(' ')[0]
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + uniprot_ns + uni_term + " ;\n"
                if 'Oryzabase' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + oryzabase_ns + records['attributes']['Oryzabase'] + " ;\n"
                if 'Oryzabase Gene Name Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['Oryzabase Gene Name Synonym(s)'].split(","):
                        syn_term = re.sub(r'"|^\s', '', syn_term)
                        os_japonica_buffer += "\t" + skos_ns + "altLabel" + "\t" + "\"%s" % (syn_term) + "\" ;\n"
                if 'Oryzabase Gene Symbol Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['Oryzabase Gene Symbol Synonym(s)'].split(","):
                        syn_term = re.sub(r'"|^\s', '', syn_term)
                        os_japonica_buffer += "\t" + skos_ns + "altSymbol" + "\t" + "\"%s" % (syn_term) + "\" ;\n"
                        #os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['attributes']['Oryzabase Gene Symbol Synonym(s)']) + " ;\n"
                if 'RAP-DB Gene Name Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['RAP-DB Gene Name Synonym(s)'].split(","):
                        syn_term = re.sub('"|\s', '', syn_term)
                        os_japonica_buffer += "\t" + skos_ns + "altLabel" + "\t" + "\"%s" % (syn_term) + "\" ;\n"
                if 'RAP-DB Gene Symbol Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['RAP-DB Gene Symbol Synonym(s)'].split(","):
                        syn_term = re.sub(r'"|^\s', '', syn_term)
                        os_japonica_buffer += "\t" + skos_ns + "altSymbol" + "\t" + "\"%s" % (syn_term) + "\" ;\n"
                        #os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['attributes']['RAP-DB Gene Symbol Synonym(s)']) + " ;\n"
                if 'Transcript_evidence' in records['attributes']:
                    for gene_id in records['attributes']['Transcript_evidence'].split(","):
                        if not (gene_id == " "):
                            if (gene_id[-1] == '.'):
                                gene_id = re.sub('.$', '', gene_id)
                            if ncbi_pattern.match(gene_id):
                                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ncbi_gene_ns + gene_id.split(" ")[0] + " ;\n"
                if 'NIAS_FLcDNA' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "evidence" + "\t" + "\"NIAS_FLcDNA:%s" % (records['attributes']['NIAS_FLcDNA']) + "\" ;\n"
                if 'TENOR' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + tenor_ns + records['attributes']['TENOR']  + " ;\n"
                if 'Expression' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"Expression:%s" % (records['attributes']['Expression']) + "\" ;\n"
                if 'B5toI1' in records['attributes']:
                    os_japonica_buffer += "\t" + skos_ns + "note" + "\t" + "\"%s" % (records['attributes']['B5toI1']) + "\" ;\n"
                if 'Manual Curation' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "hasCuration" + "\t" + "\"%s" % (records['attributes']['Manual Curation']) + "\" ;\n"
                if 'KEGG' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"dosa:%s" % (records['attributes']['KEGG']) + "\" ;\n"
                if 'Comment' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"%s" % (records['attributes']['Comment']) + "\" ;\n"

                os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' +\
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "polypeptide":
                os_japonica_buffer = ''
                os_japonica_buffer += protein_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000104" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "derives_from" + \
                                      "\t\t" + rapdb_mrna_ns + records['attributes']['Derives_from'] + " ; \n"
                os_japonica_buffer += "\t" + faldo_ns + "location" + \
                                      "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" \
                                      + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "CDS":
                os_japonica_buffer = ''
                number_cds = getCDSNumber(records['attributes']['Parent'],records['start'], number_cds)
                os_japonica_buffer += res_ns + records['attributes']['Parent'] + "#CDS" + str(number_cds) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "CDS" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000316" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "partOf" + "\t\t" + rapdb_mrna_ns + records['attributes']['Parent'] + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + "" + chromosome_ns + re.sub('Os', '', records['seqid']) + " .\n"

                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "exon":
                os_japonica_buffer = ''
                number_exon += 1
                os_japonica_buffer += res_ns + records['attributes']['Parent'] + "#exon" + str(number_exon) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Exon" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "partOf" + "\t\t" + rapdb_mrna_ns + records['attributes']['Parent'] + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "three_prime_UTR":
                os_japonica_buffer = ''
                #number_three_prime_UTR += 1
                # only one 3_prime_UTR
                #os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#three_prime_UTR_" + str(number_three_prime_UTR) + "\n"
                os_japonica_buffer += res_ns + records['attributes']['Parent'] + "#three_prime_UTR" + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "three_prime_UTR" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000205" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "partOf" + "\t\t" + rapdb_mrna_ns + records['attributes']['Parent'] + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + " " + chromosome_ns + re.sub('Os', '', records['seqid']) + " .\n"
                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "five_prime_UTR":
                os_japonica_buffer = ''
                #number_five_prime_UTR += 1
                #os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#five_prime_UTR_" + str(number_five_prime_UTR) + "\n"
                os_japonica_buffer += res_ns + records['attributes']['Parent'] + "#five_prime_UTR" + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "five_prime_UTR" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000204" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "partOf" + "\t\t" + rapdb_mrna_ns + records['attributes']['Parent'] + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + " " + chromosome_ns + re.sub('Os', '', records['seqid']) + " .\n"

                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

                # Region
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']]['number'] + ':' + \
                                      str(records['start']) + '-' + str(records['end']) + ":" + strand + "  \n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + "\";\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['start']) + ":" + strand + "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ":" + str(
                    records['end']) + ":" + strand + "  .\n\n"

                # Position 1
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['start']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['start']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

                # Position 2
                os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + chromosome_dict[records['seqid']][
                    'number'] + ":" + str(records['end']) + ":" + strand
                os_japonica_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
                os_japonica_buffer += "  ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "position" + "\t" + str(records['end']) + " ;\n"
                os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                    ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)
    strand = "1"
    position = "ForwardStrandPosition"
    print(line_number)

def annotation2RDF(path, path_output):
    print("*************Parsing Rapdb gene annotation data ***********\n")
    array = pd.read_csv(path, sep='\t', delimiter=None, dtype='str', skip_blank_lines=True)
    array.fillna('', inplace=True)

    row_count = sum(1 for row in array)
    print("Number of genes: %s\n" % (str(row_count)))
    print("Rapdb gene annotation data has been parsed!\n")
    print("*************************************\n\n")

    print("************* OryzaBase RDF conversion begins***********\n")

    ttl_handle = open(path_output, "w")

    ttl_handle.write(str(getRDFHeaders()))
    for records in array.as_matrix(columns=None):
        taxon_id = "39947"

    #    pp.pprint(orygene_ds)

        os_japonica_buffer = ''
        os_japonica_buffer += "<"+ base_resource_uri +  "transcript/" + records[0] + ">\n"
        os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
        os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records[0] + "\" ;\n"
        os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        os_japonica_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + records[0] + "\" ;\n"
        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ensembl_transcript_ns + records[0] + ";\n"
        os_japonica_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t\t" + res_ns + records[1] + " ;\n"
        os_japonica_buffer += "\t" + dc_ns + "description" + "\t" + "\"%s" % (re.sub('\"|\'', '', str(records[2]))) + "\" ;\n"
        if records[3]:
            records[3] = re.sub('\"|\'', '', str(records[3]))
            term_list = records[3].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (uni_term) + " ;\n"
        if records[4]:
            records[4] = re.sub('\"|\'', '', str(records[4]))
            term_list = records[4].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (re.sub('\"|\'', '', str(uni_term))) + " ;\n"
        if records[5]:
            records[5] = re.sub('\"|\'', '', str(records[5]))
            term_list = records[5].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "prefSymbol" + "\t" + '"%s"' % (uni_term) + " ;\n"
        if records[6]:
            records[6] = re.sub('\"|\'', '', str(records[6]))
            term_list = records[6].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "prefLabel" + "\t" + '"%s"' % (re.sub('\"|\'', '', str(uni_term))) + " ;\n"
        if records[7]:
            records[7] = re.sub('\"|\'', '', str(records[7]))
            term_list = records[7].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (uni_term) + " ;\n"
        if records[8]:
            records[8] = re.sub('\"|\'', '', str(records[8]))
            term_list = records[8].split(',')
            for uni_term in term_list:
                uni_term = re.sub('^\s+|\s+$', '', uni_term)
                os_japonica_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (re.sub('\"|\'', '', str(uni_term))) + " ;\n"
        if records[9]: # GO
            for go_term in re.findall(r'GO:[0-9]{7}', records[9]):
                os_japonica_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + re.sub(':', '_',
                                                                                                go_term) + " ;\n"
        if records[10]: # interpro
            for ipr_term in re.findall(r'IPR[0-9]{6}', records[10]):
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + interpro_ns + ipr_term + " ;\n"
        if records[11]: # transcript evidence
            if "(DDBJ, Best hit)" in records[11]:
                term = records[11].split(' ')[0]
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ena_embl_ns + re.sub('\s+', '', term) + " ;\n"
        if records[12]: # ORF evidence
            if '(UniProt)' in records[12]:
                uni_term = re.split('\s|\(', records[12])[0]
                uni_term = re.sub('\s+', '', uni_term)
                uni_term = re.sub(',', '', uni_term)
                if re.match(prot_pattern, uni_term):
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + uniprot_ns + uni_term + " ;\n"
            elif records[12].split(','):
                term_list = records[12].split(',')
                for uni_term in term_list:
                    uni_term =  re.split('\s|\(',records[12] )[0]
                    uni_term = re.sub('\s+', '', uni_term)
                    uni_term = re.sub(',', '', uni_term)
                    if re.match(prot_pattern, uni_term):
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + uniprot_ns + uni_term + " ;\n"
        if records[13]: # curration date
            os_japonica_buffer += "\t" + base_vocab_ns + "curationDate" + "\t\t" + "\"" + records[13] + "\" ;\n"
        if records[14]: # PMID
            if pubmed_pattern.match(records[14]):
                os_japonica_buffer += "\t" + dc_ns + "references" + "\t\t" + pubmed_ns + records[14] + " ;\n"
        #if records[15]: # FLcdna
        if records[16]: # oryzabase ID
            os_japonica_buffer += "\t" + owl_ns + "sameAs" + "\t\t" + oryzabase_ns + records[16] + " ;\n"

        os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
        ttl_handle.write(os_japonica_buffer)

        print(os_japonica_buffer)

    ttl_handle.close()
    print("************* Rapdb Annotation RDF completed ************!\n\n")
pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
path = '/Users/plarmande/workspace2015/datasets/head.gff'
path_output = '/Users/plarmande/workspace2015/datasets/head.ttl' # The output

#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
'''
parse the GFF from predicted and evidence from RAPDB
'''
ds = parseGFF3(path)  # call the function in another package
pp.pprint(ds)    # For to see in teminal the parsing
'''
print the RDF from GFF file parsed
'''
rapdbModeleRDF(ds, path_output)
#annotation2RDF(path,path_output)


# guideline to build resource uri
# domain name: http://www.southgreen.fr/agrold/
# resource: resource
# database : rapdb
# release-version : 5 (integer)
# species_name: oryza_sativa
# annotation_project : IRGSP-1.0
# entity_type : chromosome
# example : <http://www.southgreen.fr/agrold/resource/rapdb/5/oryza_sativa/IRGSP-1.0/chromosome/IRGSP-1.0:1:1-43270923:1
