import sys
print (sys.path)
from globalVars import *
from globalVars import base_vocab_ns
from gffParser import *
import pprint
import re
import os

'''
Created on May, 2017
The rapdbParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling rapdb data

TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: larmande
'''
__author__  = "larmande"


def getStrandValue(strandVar):
    if strandVar == "-":
        strandVar = "-1"
        positionVar = "ReverseStrandPosition"
    else:
        strandVar = "1"
        positionVar = "ForwardStrandPosition"
    return (strandVar,positionVar)

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



    pubmed_pattern = re.compile(r'^\d+$')
    ncbi_pattern = re.compile(r'^[A-Z]{2}\d{6}$')

    chromosome_size = [43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856]
    # chromosome_size.reverse()
    # The first wrinting in the file is the prefix


    print ("*************RapDB RDF conversion begins***********\n")
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    rdf_writer.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    rdf_writer.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
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
    rdf_writer.write(pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n")
    # Ajout du prefix pour la realese des donnees
    rdf_writer.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n\n")
    rdf_writer.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n\n")

# In here we buil the modele and writer in file with ttl format

  # Species triple


    os_japonica_buffer = ''
    os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + rdfs_ns + "subClassOf" + "\t\t" + sio_ns + "SIO_000253" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "OBI_0100026" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + skos_ns + "prefLabel" + "\t\t" + "Oryza sativa Japonica Group" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + skos_ns + "altLabel" + "\t\t" + "Japanese rice" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + dc_ns + "identifier" + "\t\t" + "39947" + " .\n\n"
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
            chromosome_dict[records['seqid']] = {}
            chromosome_dict[records['seqid']] = { 'uri': "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1", 'seqid': records['seqid'], 'number' : str(ch_number), 'nucleotide' : str(chromosome), 'assembly':'IRGSP-1.0'}
            chromosome_list[records['seqid']]= "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1"
            os_japonica_buffer += chromosome_ns + "IRGSP-1.0:" + str(ch_number) + ":1-" + str(chromosome) + ":1" + "\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Chromosome" + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + "Oryza sativa Japonica Group chromosome:" + "IRGSP-1.0:" + str(
                ch_number) + ":1-" + str(chromosome) + ":1" + " (IRGSP-1.0)" + "\" ;\n"
            os_japonica_buffer += "\t" + dc_ns + "identifier " + "\t" + " \"" + "IRGSP-1.0:" + str(
                ch_number) + ":1-" + str(chromosome) + ":1" + "\" ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "genomeAssembly " + "\t" + " \"" + "IRGSP-1.0" + "\" ;\n"

            print(os_japonica_buffer)
            # os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
            rdf_writer.write(os_japonica_buffer)

        #RAPDB.gff3

        if records['source']:
            if records['type'] == "gene":
                # < http: // rdf.ebi.ac.uk/resource/ ensembl / Os01g0963000 > rdfs: seeAlso < http: // identifiers.org / ensembl / Os01g0963000 >.
                #< http: // identifiers.org/ensembl / Os01g0963000 > rdf: type identifiers: ensembl.
                # < http: // identifiers.org/ensembl/Os01g0963000> sio: SIO_000671[a ident_type: ensembl; sio: SIO_000300 "Os01g0963000"].
               #  < http: // rdf.ebi.ac.uk/resource/ensembl/Os01g0963000> faldo: location < http: // rdf.ebi.ac.uk / resource / ensembl / 41 / oryza_sativa / IRGSP - 1.0 / 1: 42441660 - 42442936:1 >.
               #  strand = "1"
               #  position = "ForwardStrandPosition"
               #  if records['strand'] == "-":
               #      strand = "-1"
               #      position = "ReverseStrandPosition"
                (strand,position) = getStrandValue(records['strand'])

                os_japonica_buffer = ''

                os_japonica_buffer += rapdb_gene_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" +  ensembl_gene_ns + records['attributes']['ID'] + " ;\n"
                os_japonica_buffer += "\t" + owl_ns + "sameAs" + "\t" + ensembl_gene_ns + records['attributes']['ID'] + " ;\n"

                # rapdb..ID  skos:closeMatch ensembl:id  ## important ##
                os_japonica_buffer += "\t" + base_vocab_ns + "source_project" + "\t" + " \"" + records['source'] + "\" ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Gene" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_biotype" + "\t" + "\"protein_coding\" ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                os_japonica_buffer += "\t" + dc_ns + "identifier" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"

                os_japonica_buffer += "\t" + dc_ns + "description" + "\t" + " \"" + records['attributes']['Note'] + "\" ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
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

                ## <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701-42442930:1> rdfs:label "chromosome 1:42441701-42442930:1" .
                ## <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701-42442930:1> rdf:type faldo:Region .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701-42442930:1> faldo:begin <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701:1> .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701-42442930:1> faldo:end <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42442930:1> .
                # question #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701-42442930:1> faldo:reference <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1> .

                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701:1> rdf:type faldo:ExactPosition .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701:1> rdf:type faldo:ForwardStrandPosition .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701:1> faldo:position 42441701 .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42441701:1> faldo:reference <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1> .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42442930:1> rdf:type faldo:ExactPosition .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42442930:1> rdf:type faldo:ForwardStrandPosition .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42442930:1> faldo:position 42442930 .
                #<http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1:42442930:1> faldo:reference <http://rdf.ebi.ac.uk/resource/ensembl/41/oryza_sativa/IRGSP-1.0/1> .

                # < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01 > sio: SIO_000974 < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01  # Exon_1> .
                # < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01  # Exon_1> rdf:type sio:SIO_001261 .
                # < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01  # Exon_1> sio:SIO_000628 <http://rdf.ebi.ac.uk/resource/ensembl.exon/Os01t0963000-01.exon1> .
                # < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01  # Exon_1> sio:SIO_000300 1 .
                # < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 > obo: RO_0002162 taxon: 39947.
                # < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 > dc: identifier "Os01t0963000-01".
                # < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 > rdfs: seeAlso < http: // identifiers.org / ensembl / Os01t0963000 - 01 >.
                # < http: // identifiers.org / ensembl / Os01t0963000 - 01 > rdf: type identifiers: ensembl.
                # < http: // identifiers.org / ensembl / Os01t0963000 - 01 > sio: SIO_000671[a ident_type: ensembl; sio: SIO_000300 "Os01t0963000-01"].
                # < http: // rdf.ebi.ac.uk / resource / ensembl.transcript / Os01t0963000 - 01 > obo: SO_translates_to < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 >.
                # < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 > rdf: type term: protein.
                # < http: // rdf.ebi.ac.uk / resource / ensembl.protein / Os01t0963000 - 01 > rdfs: seeAlso gene3d: 1.10.420.10.
                # ...
                # describing the ressource gene < http: // rdf.ebi.ac.uk / resource / ensembl / Os01g0963000 > sio: SIO_000558 ensembl: Os01g0962900.
                # <http://rdf.ebi.ac.uk/resource/ensembl.transcript/Os01t0963000-01> obo:SO_transcribed_from <http://rdf.ebi.ac.uk/resource/ensembl/Os01g0963000> .
                # <http://rdf.ebi.ac.uk/resource/ensembl.protein/Os01t0963000-01> rdf:type term:protein .
                # <http://rdf.ebi.ac.uk/resource/ensembl.transcript/Os01t0963000-01> obo:SO_translates_to <http://rdf.ebi.ac.uk/resource/ensembl.protein/Os01t0963000-01> .
                # <http://rdf.ebi.ac.uk/resource/ensembl/Os01g0963000> sio:SIO_000558 ensembl:Os01g0962900 . # is orthologuos to
                print(os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)


            if records['type'] == "mRNA":

                (strand,position) = getStrandValue(records['strand'])
                os_japonica_buffer = ''
                os_japonica_buffer += rapdb_mrna_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "mRNA" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "source_project" + "\t" + " \"" + records['source'] + "\" ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_biotype" + "\t" + "\"protein_coding\" ;\n"
                os_japonica_buffer += "\t" + dc_ns + "identifier" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ensembl_transcript_ns + records['attributes']['Name'] + ";\n"

                if 'Note' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "comment" + "\t" + '"%s"' % (records['attributes']['Note']) + " ;\n"

                # os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "develops_from" + "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"
                os_japonica_buffer += "\t" + obo_ns +  "SO_transcribed_from" +  "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"
                # <http://rdf.ebi.ac.uk/resource/ensembl.transcript/Os09t0372700-01> obo:SO_translates_to <http://rdf.ebi.ac.uk/resource/ensembl.protein/Os09t0372700-01> .
                ## os_japonica_buffer += "\t" + obo_ns + "SO_translates_to" + "\t\t" + rapdb_gene_ns + records['attributes']['Locus_id'] + " ;\n"

                os_japonica_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + "IRGSP-1.0:" + \
                                      chromosome_dict[records['seqid']]['number'] + ':' + str(
                    records['start']) + '-' + str(records['end']) + ":" + strand + " .\n\n"

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

                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + "" + chromosome_ns + re.sub('Os', '', records['seqid']) + " ;\n"

                if 'GO' in records['attributes']:
                        for go_term in re.findall(r'GO:[0-9]{7}',records['attributes']['GO']):
                            os_japonica_buffer += "\t" + base_vocab_ns + "go_term" + "\t" + obo_ns + re.sub(':', '_', go_term) + " ;\n"
                            #os_japonica_buffer += "\t" + base_vocab_ns + "comment" + "\t" + '"%s"' % (records['attributes']['GO']) + " ;\n"
                if 'InterPro' in records['attributes']:
                    for ipr_term in re.findall(r'IPR[0-9]{6}', records['attributes']['InterPro']):
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" +"\t" + interpro_ns + ipr_term + " ;\n"
                    os_japonica_buffer += "\t" + base_vocab_ns + "comment" + "\t" + '"%s"' % (records['attributes']['InterPro']) + " ;\n"
                if 'CGSNL Gene Name' in records['attributes']:
                    syn_term = re.sub('"', '', records['attributes']['CGSNL Gene Name'])
                    os_japonica_buffer += "\t" + base_vocab_ns + "has_synonym" + "\t" + '"%s"' % (syn_term) + " ;\n"
                if 'CGSNL Gene Symbol' in records['attributes']:
                    sym_term = re.sub('"', '', records['attributes']['CGSNL Gene Symbol'])
                    os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (sym_term) + " ;\n"
                if 'Literature_PMID' in records['attributes']:
                    if pubmed_pattern.match(records['attributes']['Literature_PMID']):
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + pubmed_ns + records['attributes']['Literature_PMID'] + " ;\n"
                if 'ORF_evidence' in records['attributes']:
                    if '(UniProt)' in records['attributes']['ORF_evidence']:
                        uni_term = records['attributes']['ORF_evidence'].split(' ')[0]
                        os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + uniprot_ns + uni_term + " ;\n"
                if 'Oryzabase' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + oryzabase_ns + records['attributes']['Oryzabase'] + " ;\n"
                if 'Oryzabase Gene Name Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['Oryzabase Gene Name Synonym(s)'].split(","):
                        syn_term = re.sub('"', '', syn_term)
                        os_japonica_buffer += "\t" + base_vocab_ns + "has_synonym" + "\t" + '"%s"' % (syn_term) + " ;\n"
                if 'Oryzabase Gene Symbol Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['Oryzabase Gene Symbol Synonym(s)'].split(","):
                        syn_term = re.sub('"', '', syn_term)
                        os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (syn_term) + " ;\n"
                        #os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['attributes']['Oryzabase Gene Symbol Synonym(s)']) + " ;\n"
                if 'RAP-DB Gene Name Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['RAP-DB Gene Name Synonym(s)'].split(","):
                        syn_term = re.sub('"', '', syn_term)
                        os_japonica_buffer += "\t" + base_vocab_ns + "has_synonym" + "\t" + '"%s"' % (syn_term) + " ;\n"
                if 'RAP-DB Gene Symbol Synonym(s)' in records['attributes']:
                    for syn_term in records['attributes']['RAP-DB Gene Symbol Synonym(s)'].split(","):
                        syn_term = re.sub('"', '', syn_term)
                        os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (syn_term) + " ;\n"
                        #os_japonica_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['attributes']['RAP-DB Gene Symbol Synonym(s)']) + " ;\n"
                if 'Transcript_evidence' in records['attributes']:
                    for gene_id in records['attributes']['Transcript_evidence'].split(","):
                        if not (gene_id == " "):
                            if (gene_id[-1] == '.'):
                                gene_id = re.sub('.$', '', gene_id)
                            if ncbi_pattern.match(gene_id):
                                os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ncbi_gene_ns + gene_id.split(" ")[0] + " ;\n"
                if 'NIAS_FLcDNA' in records['attributes']:
                    os_japonica_buffer += "\t" + base_vocab_ns + "evidence" + "\t" + '"NIAS_FLcDNA: %s"' % (records['attributes']['NIAS_FLcDNA']) + " ;\n"
                if 'TENOR' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + tenor_ns + records['attributes']['TENOR']  + " ;\n"
                if 'Expression' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + '"Expression: %s"' % (records['attributes']['Expression']) + " ;\n"
                if 'B5toI1' in records['attributes']:
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + '"B5toI1: %s"' % (records['attributes']['B5toI1']) + " ;\n"
                print(os_japonica_buffer)
                os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
                rdf_writer.write(os_japonica_buffer)

            if records['type'] == "polypeptide":
                os_japonica_buffer = ''
                os_japonica_buffer += protein_ns + records['attributes']['ID'] + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Protein" + " ;\n"
                #os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000104" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "derives_from" + "\t\t" + mRNA_ns + records['attributes']['Derives_from'] + " ; \n"
                # os_japonica_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t\t" + "" + chromosome_ns + re.sub('Os', '', records['seqid']) + " .\n"
                # os_japonica_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "IRGSP-1.0:" + str(
                #     ch_number) + ":1-" + str(chromosome) + ":1" + " .\n\n"

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

            if records['type'] == "CDS":
                os_japonica_buffer = ''
                number_cds += 1
                os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#CDS" + str(number_cds) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "CDS" + " ;\n"
                #os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000316" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "part_of" + "\t\t" + mRNA_ns + records['attributes']['Parent'] + " ;\n"
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
                os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#exon" + str(number_exon) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Exon" + " ;\n"
                #os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000147" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "part_of" + "\t\t" + mRNA_ns + records['attributes']['Parent'] + " ;\n"
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

            if records['type'] == "three_prime_UTR":
                os_japonica_buffer = ''
                number_three_prime_UTR += 1
                os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#three_prime_UTR_" + str(number_three_prime_UTR) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Threee_prime_UTR" + " ;\n"
                #os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000205" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "part_of" + "\t\t" + mRNA_ns + records['attributes']['Parent'] + " ;\n"
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
                number_five_prime_UTR += 1
                os_japonica_buffer += OrygenesDB_ns + records['attributes']['Parent'] + "#five_prime_UTR_" + str(number_five_prime_UTR) + "\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + res_ns + "Five_prime_UTR" + " ;\n"
                #os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                #os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000204" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
                os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + " \"" + str(records['start']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + " \"" + str(records['end']) + "\"^^xsd:integer ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "part_of" + "\t\t" + mRNA_ns + records['attributes']['Parent'] + " ;\n"
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


pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
path = '/Users/plarmande/Downloads/IRGSP-1.0_representative_12-18/test.all.gff'
path_output = '/Users/plarmande/Downloads/IRGSP-1.0_representative_12-18/Oryza_sativa_Japonica.ttl' # The output

#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
ds = parseGFF3(path)   # The parsing file with the tropGeneParser()
# pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

rapdbModeleRDF(ds, path_output)

# guideline to build resource uri
# domain name: http://www.southgreen.fr/agrold/
# resource: resource
# database : rapdb
# release-version : 5 (integer)
# species_name: oryza_sativa
# annotation_project : IRGSP-1.0
# entity_type : chromosome
# example : <http://www.southgreen.fr/agrold/resource/rapdb/5/oryza_sativa/IRGSP-1.0/chromosome/IRGSP-1.0:1:1-43270923:1