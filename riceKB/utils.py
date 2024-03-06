import sys
print(sys.path)
from globalVars import *
from globalVars import base_vocab_ns
from gffParser import *
import rdflib
from rdflib.graph import Graph
import pprint
import re
import os

'''
Created on Dec, 2019
The utils module is created as part of the AgroLD project.


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

def getChromosomeNumber(chrString):
    '''
    returns the formated chromosome number without Chr string or others
    '''
    if 'Chr' in chrString:
        chromosome_nb = re.sub('Chr', '', chrString)
        chromosome_nb = re.sub('^0', '', chromosome_nb)
    else:
        chromosome_nb = re.sub('^0', '',chrString)
    return chromosome_nb

def getFaldoRegion(taxon_id,ssp, seqid,start,end,strand):
    '''
    input taxon_if, sub-species name, Chromomose ID, Start, end, strand
    returns a formated block of Faldo region, positions
    '''
    (strand, position) = getStrandValue(strand)
    strand = str(strand)
    genome_buffer =''
    # Region

    genome_buffer += "<" + chromosome_uri + taxon_id + "/" + ssp \
                     +  seqid + ":" + \
                     str(start) + "-" + str(end) + ":" \
                     + strand + "> \n"
    genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + taxon_id +"/"+ ssp \
                     +  seqid + ":" + \
                     str(start) + "-" + str(end) + ":" \
                     + strand + "\";\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
    genome_buffer += "\t" + faldo_ns + "begin" + "\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     + seqid + ":" + \
                     str(start) + ":" + strand + ">;\n"
    genome_buffer += "\t" + faldo_ns + "end" + "\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     + seqid + ":" + \
                     str(end) + ":" + strand + ">  .\n\n"

    # Position 1
    genome_buffer += "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     +  seqid + ":" + \
                     str(start) + ":" + strand + ">\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
    genome_buffer += "  ;\n"
    genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(start) + " ;\n"
    genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     +  seqid + "> .\n\n"

    # Position 2
    genome_buffer += "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     +  seqid + ":" + \
                     str(end) + ":" + strand + "> \n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
    genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
    genome_buffer += " ;\n"
    genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(end) + " ;\n"
    genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                     +  seqid + "> .\n\n"
    return genome_buffer

# def getFaldoRegion(taxon_id,seqid,start,end,strand):
#     (strand, position) = getStrandValue(strand)
#     strand = str(strand)
#     genome_buffer =''
#     # Region
#
#     genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(start) + "-" + str(end) + ":" \
#                      + strand + "> \n"
#     genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(start) + "-" + str(end) + ":" \
#                      + strand + "\";\n"
#     genome_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
#     genome_buffer += "\t" + faldo_ns + "begin" + "\t" + "<" + chromosome_uri + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(start) + ":" + strand + ">;\n"
#     genome_buffer += "\t" + faldo_ns + "end" + "\t" + "<" + chromosome_uri + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(end) + ":" + strand + ">  .\n\n"
#
#     # Position 1
#     genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(start) + ":" + strand + ">\n"
#     genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
#     genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
#     genome_buffer += "  ;\n"
#     genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(start) + " ;\n"
#     genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + \
#                      "/" + re.sub('Os|Chr', '', seqid) + "> .\n\n"
#
#     # Position 2
#     genome_buffer += "<" + chromosome_uri + taxon_id + "/" \
#                      + re.sub('Os|Chr', '', seqid) + ":" + \
#                      str(start) + ":" + strand + "> \n"
#     genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
#     genome_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
#     genome_buffer += " ;\n"
#     genome_buffer += "\t" + faldo_ns + "position" + "\t" + str(end) + " ;\n"
#     genome_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + \
#                      "/" + re.sub('Os|Chr', '', seqid) + "> .\n\n"
#     return genome_buffer

def getRDFHeaders():
    headersBuffer = ''
    headersBuffer += base + "\t" + "<" + base_uri + "> .\n"
    headersBuffer += pr + "\t" + rdf_ns + "<" + rdf + "> .\n"
    headersBuffer += pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n"
    headersBuffer += pr + "\t" + xsd_ns + "<" + xsd + "> .\n"
    headersBuffer += pr + "\t" + owl_ns + "<" + owl_uri + "> .\n"
    headersBuffer += pr + "\t" + dc_ns + "<" + dc_uri + "> .\n"
    headersBuffer += pr + "\t" + dcterms_ns + "<" + dcterms_uri + "> .\n"
    headersBuffer += pr + "\t" + skos_ns + "<" + skos + "> .\n"
    headersBuffer += pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n"
    headersBuffer += pr + "\t" + obo_ns + "<" + obo_uri + "> .\n"
    headersBuffer += pr + "\t" + sio_ns + "<" + sio_uri + "> .\n"
    headersBuffer += pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n"
    headersBuffer += pr + "\t" + interpro_ns + "<" + interpro_uri + "> .\n"
    headersBuffer += pr + "\t" + tenor_ns + "<" + tenor_uri + "> .\n"
    headersBuffer += pr + "\t" + ncbi_gene_ns + "<" + ncbi_gene_uri + "> .\n"
    headersBuffer += pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n"
    headersBuffer += pr + "\t" + faldo_ns + "<" + faldo + "> .\n"
    headersBuffer += pr + "\t" + res_ns + "<" + resource + "> .\n"
    headersBuffer += pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n"
    headersBuffer += pr + "\t" + ena_embl_ns + "<" + ena_embl_uri + "> .\n"
    headersBuffer += pr + "\t" + pubmed_ns + "<" + pubmed_uri + "> .\n"
    headersBuffer += pr + "\t" + pubmed_rdf_ns + "<" + pubmed_rdf_uri + "> .\n"
    headersBuffer += pr + "\t" + ensembl_gene_ns + "<" + ensembl_gene_uri + "> .\n"
    headersBuffer += pr + "\t" + ensembl_transcript_ns + "<" + ensembl_transcript_uri + "> .\n"
    headersBuffer += pr + "\t" + ensembl_protein_ns + "<" + ensembl_protein_uri + "> .\n"
    headersBuffer += pr + "\t" + uniprot_ns + "<" + uniprot_uri + "> .\n"
    headersBuffer += pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n"
    headersBuffer += pr + "\t" + rapdb_gene_ns + "<" + rapdb_gene_uri + "> .\n"
    headersBuffer += pr + "\t" + rapdb_mrna_ns + "<" + rapdb_mrna_uri + "> .\n"
    headersBuffer += pr + "\t" + gr_g_ns + "<" + gramene_gene + "> .\n"
    headersBuffer += pr + "\t" + up_ns + "<" + uniprot + "> .\n"
    headersBuffer += pr + "\t" + bibo_ns + "<" + bibo_uri + "> .\n"
    headersBuffer += pr + "\t" + prism_ns + "<" + prism_uri + "> .\n"
    # Ajout du prefix pour la realese des donnees
    headersBuffer += pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n\n"
    return headersBuffer

def RDF_validation(ttl_buffer,ttl_handle,oryid,output_dir):

    try:
        temp_file = os.path.join(output_dir, 'temp_graph.ttl')
        #temp_file = '/Users/pierre/Downloads/tmp/temp_graph.ttl'
        try_handle = open(temp_file, "w")
        try_handle.write(str(getRDFHeaders()))
        try_handle.write(ttl_buffer)
        try_handle.close()
        #output_file = '/Users/plarmande/Downloads/oryzabase.nt'
        g = Graph()
        g.parse(temp_file, format="turtle")
        # output = open(output_file, "w")
        # w = Graph()
        #g.serialize(destination="file:/Users/plarmande/Downloads/oryzabase_test.ttl", format='xml')
        ttl_handle.write(ttl_buffer)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        temp = os.path.join(output_dir, 'temp_graph' + oryid + '.ttl')
        #temp = '/Users/pierre/Downloads/tmp/temp_graph'+ oryid +'.ttl'
        handle = open(temp, "w")
        handle.write(str(getRDFHeaders()))
        handle.write(ttl_buffer)
        pass
