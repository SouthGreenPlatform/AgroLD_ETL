from __future__ import print_function
import sys
#print sys.path
#from riceKB.globalVars import *
from riceKB.globalVars import *
import pprint
import re
import os
# TODO modifier dans la base les predicats develops_from et has_trait avec developsFrom hasTrait
'''
Created on Sept, 2019
The ensmbl.rdf module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Esembl RDF datasets to normalize with AgroLD

TODO:
    1) Add documentation
    2) better Error handling
    
Questions: 
- which URI to define for resources such as genes, transcript and proteins ? \
for now its all after http://www.southgreen.fr/agrold/resource/ 

@author: larmande
'''
# <http://biohackathon.org/resource/faldo#begin>
# <http://biohackathon.org/resource/faldo#end>
# <http://biohackathon.org/resource/faldo#location>
# <http://biohackathon.org/resource/faldo#position>
# <http://biohackathon.org/resource/faldo#reference>
# <http://purl.obolibrary.org/obo/RO_0002162> ## in taxon
# <http://purl.obolibrary.org/obo/SO_has_part>
# <http://purl.obolibrary.org/obo/SO_transcribed_from>
# <http://purl.obolibrary.org/obo/SO_translates_to>
# <http://rdf.ebi.ac.uk/terms/ensembl/inEnsemblAssembly>
# <http://rdf.ebi.ac.uk/terms/ensembl/inEnsemblSchemaNumber>
# dc:description
# dc:identifier
# void:subset
# <http://semanticscience.org/resource/SIO_000300>  ## has value  : subject(string) > gene id
# <http://semanticscience.org/resource/SIO_000558>  ## 	is orthologous to  : geneUri > geneURI
# <http://semanticscience.org/resource/SIO_000628>  ## refers to : exon > transcript
# <http://semanticscience.org/resource/SIO_000671>  ## has identifier  : geneURI > subject(string) <http://identifiers.org/ensembl/Os04t0308600-00>	t1757965
# <http://semanticscience.org/resource/SIO_000974>  ##	has ordered part  :  transcript > exon
# rdfs:seeAlso
# rdfs:label
# rdfs:subClassOf
# owl:sameAs
# rdf:type
# skos:altLabel
# skos:altlabel ## Attention a corriger
# skos:prefLabel
# <http://rdf.ebi.ac.uk/terms/ensembl/ANNOTATED> ## ?geneURI > other URI such as Uniprot, PFAM etc. istance of <http://rdf.ebi.ac.uk/terms/ensembl/EnsemblDBEntry>
# <http://rdf.ebi.ac.uk/terms/ensembl/CHECKSUM> ## useless for us
# <http://rdf.ebi.ac.uk/terms/ensembl/DEPENDENT>
# <http://rdf.ebi.ac.uk/terms/ensembl/DIRECT>
# <http://rdf.ebi.ac.uk/terms/ensembl/INFERRED_FROM_TRANSCRIPT>
# <http://rdf.ebi.ac.uk/terms/ensembl/INFERRED_FROM_TRANSLATION>
# <http://rdf.ebi.ac.uk/terms/ensembl/SEQUENCE_MATCH>


__author__  = "larmande"

ROOT_DIR = ''
file_input = ''
prefixes = {}
patern1 = "^<http://www.southgreen.fr/agrold/resource/>"
def ensemblParser(files, type):
    """
    parse the ensembl rdf file
    :param files: path of the file
    write  objects corresponding to prefix and triples
    """

    if(os.path.isfile(str(files))):
        print("***************** Parsing Esembl RDF data in " + files + " ********************\n")
        ttl_handle = open(ensembl_out, "w")
        with open(files) as fp:
            for line in fp:
                if(re.match(r'^@prefix', line)):
                    #print line
                    prefix, value = re.split(':\s+', line)
                    prefix = re.sub('@prefix\s+', '', prefix)
                    value = re.sub('\.\n$', '', value)
                    prefixes[prefix]=value
                    if prefix == 'dataset':
                        value = '<http://www.southgreen.fr/agrold/dataset/>'
                    if prefix == 'ensembl':
                        value = '<http://www.southgreen.fr/agrold/resource/>'
                    if prefix == 'ensembl_variant':
                            value = '<http://www.southgreen.fr/agrold/resource/variant/>'
                    if prefix == 'ensemblvariation':
                        value = '<http://www.southgreen.fr/agrold/vocabulary/variation/>'
                    if prefix == 'exon':
                        value = '<http://www.southgreen.fr/agrold/resource/exon/>'
                    if prefix == 'protein':
                        value = '<http://www.southgreen.fr/agrold/resource/protein/>'
                    if prefix == 'term':
                        value = '<http://www.southgreen.fr/agrold/vocabulary/>'
                    if prefix == 'transcript':
                        value = '<http://www.southgreen.fr/agrold/resource/transcript/>'

                    ttl_handle.write("@prefix " + prefix + ": " + value + ".\n")
                    #print("@prefix " + prefix + ": " + value + ".\n")
                elif type == 'xref':
                    if re.findall("<http://www.southgreen.fr/agrold/resource/",line):
                        line = re.sub('http://rdf\.ebi\.ac\.uk/resource/ensembl/', 'http://www.southgreen.fr/agrold/resource/', line)
                    if re.findall("<http://www.southgreen.fr/agrold/resource/transcript/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.transcript/',
                                      '<http://www.southgreen.fr/agrold/resource/transcript/', line)
                    if re.findall("<http://www.southreen.fr/agrold/resource/protein/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.protein/',
                                      '<http://www.southreen.fr/agrold/resource/protein/', line)
                    if re.findall("<http://www.southreen.fr/agrold/resource/exon/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.exon/',
                                      '<http://www.southreen.fr/agrold/resource/exon/', line)
                    if re.findall("skos: altlabel",line):
                        line = re.sub('skos: altlabel',   'skos: altLabel', line)
                    ttl_handle.write(line)

                    # elif re.findall("rdfs:seeAlso",line):
                    #     #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                    #     ttl_handle.write(line)
                    # elif re.findall("rdfs:label",line):
                    #     #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                    #     ttl_handle.write(line)
                    # elif re.findall("dc:description",line):
                    #     #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                    #     ttl_handle.write(line)
                    # if re.findall("dc:identifier", line):
                    #         # line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                    #         ttl_handle.write(line)
                    # if re.findall("owl:sameAs",line):
                    #     #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                    #     ttl_handle.write(line)

                    # if re.findall("skos: altLabel",line):
                    # if re.findall("term:DEPENDENT", line):
                            # line = re.sub('term:DEPENDENT',   'term:dependent', line)
                            # ttl_handle.write(line)
                else :
                    if re.findall("<http://www.southgreen.fr/agrold/resource/",line):
                        line = re.sub('http://rdf\.ebi\.ac\.uk/resource/ensembl/\d*/?', 'http://www.southgreen.fr/agrold/resource/', line)
                        if re.findall('rdfs:subClassOf <http://www.southgreen.fr/agrold/resource/',line):
                            line = ''
                    # if re.findall("term:inEnsemblAssembly",line):
                    #     line = re.sub('term:inEnsemblAssembly',
                    #                   'term:inAssembly', line)
                    # if re.findall("term:inEnsemblSchemaNumber",line):
                    #     line = re.sub('term:inEnsemblSchemaNumber',
                    #                   'term:inSchemaNumber', line) # a supprimer ou modifier
                    if re.findall("<http://www.southgreen.fr/agrold/resource/transcript/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.transcript/',
                                      '<http://www.southgreen.fr/agrold/resource/transcript/', line)
                    if re.findall("<http://www.southreen.fr/agrold/resource/protein/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.protein/',
                                      '<http://www.southreen.fr/agrold/resource/protein/', line)
                    if re.findall("<http://www.southreen.fr/agrold/resource/exon/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.exon/',
                                      '<http://www.southreen.fr/agrold/resource/exon/', line)
                    if re.findall("term:protein",line):
                    if re.findall("term:Protein",line):
                        ttl_handle.write(line)
                        line = re.sub('term:protein', 'term:Protein', line)
                        line = re.sub('term:Protein', 'term:Protein', line)
                    if re.findall("term:Protein_coding",line):
                    if re.findall("term:Gene",line):
                        ttl_handle.write(line)
                        line = re.sub('term:Protein_coding', 'term:Gene', line)
                        line = re.sub('term:Gene', 'term:Gene', line)
                    if re.findall("term:protein_coding",line):
                    if re.findall("term:Protein_coding",line):
                    if re.findall("term:Gene",line):
                        ttl_handle.write(line)
                        line = re.sub('term:protein_coding', 'term:Gene', line)
                        line = re.sub('term:Protein_coding', 'term:Gene', line)
                        line = re.sub('term:Gene', 'term:Gene', line)
                    # if re.findall("term:EnsemblRegion", line):
                    #     line = re.sub('term:EnsemblRegion', 'term:Region', line)
                    if re.findall("<http://www.southgreen.fr/agrold/dataset/ensemblgenomes/", line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/dataset/ensemblgenomes/\d*/?',
                                      '<http://www.southgreen.fr/agrold/dataset/ensemblgenomes/', line)
                    # change   obo:SO_has_part >  term:has_part
                    # change   term:hasPart >  term:has_part
                    if re.findall("obo:SO_has_part", line):
                    if re.findall("term:hasPart", line):
                        ttl_handle.write(line)
                        line = re.sub('obo:SO_has_part', 'term:hasPart', line)
                        line = re.sub('term:hasPart', 'term:hasPart', line)
                    # obo:SO_transcribed_from > term:transcribed_from ou term:develops_from
                    # term:developsFrom > term:transcribed_from ou term:develops_from
                    if re.findall("obo:SO_transcribed_from", line):
                    if re.findall("term:developsFrom", line):
                        ttl_handle.write(line)
                        line = re.sub('obo:SO_transcribed_from', 'term:developsFrom', line)
                        line = re.sub('term:developsFrom', 'term:developsFrom', line)
                    # obo:SO_translates_to > term:translates_to ou term:encodes
                    # term:encodes > term:translates_to ou term:encodes
                    if re.findall("obo:SO_translates_to", line):
                    if re.findall("term:encodes", line):
                        ttl_handle.write(line)
                        line = re.sub('obo:SO_translates_to', 'term:encodes', line)
                        line = re.sub('term:encodes', 'term:encodes', line)
                    # delete ?o where <http://semanticscience.org/resource/SIO_000671>
                    # if re.findall("sio:SIO_000671 \[a ", line):
                    #     #line = re.sub('obo:SO_translates_to', 'obo:so#translates_to', line)
                    #     #line = re.sub('term:encodes', 'obo:so#translates_to', line)
                    #     line = ''
                    # taxon:39947 rdfs:subClassOf obo:OBI_0100026 . + add taxon:39947 owl:sameAs ncbiTaxon:39947
                    # delete http://www.southgreen.fr/agrold/resource/oryza_sativa/IRGSP-1.0/chromosome:IRGSP-1.0:1:1:43270923:1> rdfs:subClassOf http://www.southgreen.fr/agrold/resource/oryza_sativa/IRGSP-1.0/chromosome:IRGSP-1.0:1:1:43270923:1> .
                    # in taxon = http://purl.obolibrary.org/obo/RO_0002162 ; change taxon	taxon	ObjectProperty	SIO_000253	exact match (SIO_000253:has source) to http://purl.obolibrary.org/obo/RO_0002162
                    # <http://identifiers.org/ensembl/Os12t0534500-00-E1> rdf:type identifiers:ensembl
                    # rdfs:seeAlso panther: <http://purl.uniprot.org/panther/>
                    #print(line)
                    ttl_handle.write(line)
    else:
        print("***************** File not found ********************\n")
        print(files)


if len(sys.argv):
    filepath = sys.argv[0]
    ROOT_DIR = os.path.dirname(filepath)
    # print(ROOT_DIR)
    file_input = os.path.basename(filepath)
    # print(file_input)
else:
    ROOT_DIR = '/Users/pierre/workspace2015/datasets/'
    file_input = 'oryza_sativa.ttl'

ensembl_files =  os.path.join(ROOT_DIR + '/' + file_input)
ensembl_out = os.path.join(ROOT_DIR + '/' + 'agrold.' + file_input)


pp = pprint.PrettyPrinter(indent=4)

if (re.search(r'xrefs', ensembl_files)):
    ensemblParser(ensembl_files,type='xref') # type = 'xref' or None
    # print('xrefs OK')
else:
    ensemblParser(ensembl_files,type=None)

print("***************** Esembl RDF data ********************\n")

print("********************************************************\n\n")

# if __name__ == '__main__':
# 	link = "http://data.gramene.org/v60/genes?q=Os08g0494375"
# 	result = connectionError(link)
# 	print(result.content)