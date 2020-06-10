from __future__ import print_function
import sys
#print sys.path
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


__author__  = "larmande"

ROOT_DIR = ''
file_input = ''
prefixes = {}
patern1 = "^<http://rdf.ebi.ac.uk/resource/ensembl/>"
def ensemblParser(files, type):
    """
    parse the ensembl rdf file
    :param files: path of the file
    write  objects corresponding to prefix and triples
    """

    if(os.path.isfile(str(files))):
        print("***************** Parsing Esembl RDF data ********************\n")
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
                    print("@prefix " + prefix + ": " + value + ".\n")
                elif type == 'xref':
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl/",line):
                        line = re.sub('http://rdf\.ebi\.ac\.uk/resource/ensembl/', 'http://www.southgreen.fr/agrold/resource/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.transcript/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.transcript/',
                                      '<http://www.southgreen.fr/agrold/resource/transcript/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.protein/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.protein/',
                                      '<http://www.southreen.fr/agrold/resource/protein/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.exon/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.exon/',
                                      '<http://www.southreen.fr/agrold/resource/exon/', line)
                    if re.findall("rdfs:seeAlso",line):
                        #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                        ttl_handle.write(line)
                    if re.findall("rdfs:label",line):
                        #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                        ttl_handle.write(line)
                    if re.findall("dc:description",line):
                        #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                        ttl_handle.write(line)
                    if re.findall("dc:identifier", line):
                            # line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                            ttl_handle.write(line)
                    if re.findall("owl:sameAs",line):
                        #line = re.sub('term:inEnsemblSchemaNumber',   'term:inSchemaNumber', line)
                        ttl_handle.write(line)
                    if re.findall("term:DEPENDENT", line):
                            line = re.sub('term:DEPENDENT',   'term:dependent', line)
                            ttl_handle.write(line)
                else :
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl/",line):
                        line = re.sub('http://rdf\.ebi\.ac\.uk/resource/ensembl/\d*/?', 'http://www.southgreen.fr/agrold/resource/', line)
                        if re.findall('rdfs:subClassOf <http://www.southgreen.fr/agrold/resource/',line):
                            line = ''
                    if re.findall("term:inEnsemblAssembly",line):
                        line = re.sub('term:inEnsemblAssembly',
                                      'term:inAssembly', line)
                    if re.findall("term:inEnsemblSchemaNumber",line):
                        line = re.sub('term:inEnsemblSchemaNumber',
                                      'term:inSchemaNumber', line) # a supprimer ou modifier
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.transcript/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.transcript/',
                                      '<http://www.southgreen.fr/agrold/resource/transcript/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.protein/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.protein/',
                                      '<http://www.southreen.fr/agrold/resource/protein/', line)
                    if re.findall("<http://rdf.ebi.ac.uk/resource/ensembl.exon/",line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/resource/ensembl\.exon/',
                                      '<http://www.southreen.fr/agrold/resource/exon/', line)
                    if re.findall("term:protein",line):
                        line = re.sub('term:protein', 'term:Protein', line)
                    if re.findall("term:Protein_coding",line):
                        line = re.sub('term:Protein_coding', 'term:Gene', line)
                    if re.findall("term:protein_coding",line):
                        line = re.sub('term:protein_coding', 'term:Gene', line)
                    if re.findall("term:EnsemblRegion", line):
                        line = re.sub('term:EnsemblRegion', 'term:Region', line)
                    if re.findall("<http://rdf.ebi.ac.uk/dataset/ensemblgenomes/", line):
                        line = re.sub('<http://rdf\.ebi\.ac\.uk/dataset/ensemblgenomes/\d*/?',
                                      '<http://www.southgreen.fr/agrold/dataset/ensemblgenomes/', line)
                    # change   obo:SO_has_part >  term:has_part
                    if re.findall("obo:SO_has_part", line):
                        line = re.sub('obo:SO_has_part', 'term:hasPart', line)
                    # obo:SO_transcribed_from > term:transcribed_from ou term:develops_from
                    if re.findall("obo:SO_transcribed_from", line):
                        line = re.sub('obo:SO_transcribed_from', 'term:developsFrom', line)
                    # obo:SO_translates_to > term:translates_to ou term:encodes
                    if re.findall("obo:SO_translates_to", line):
                        line = re.sub('obo:SO_translates_to', 'term:encodes', line)
                    # delete ?o where <http://semanticscience.org/resource/SIO_000671>
                    if re.findall("sio:SIO_000671 \[a ", line):
                        #line = re.sub('obo:SO_translates_to', 'obo:so#translates_to', line)
                        line = ''
                    # taxon:39947 rdfs:subClassOf obo:OBI_0100026 . + add taxon:39947 owl:sameAs ncbiTaxon:39947
                    # delete http://www.southgreen.fr/agrold/resource/oryza_sativa/IRGSP-1.0/chromosome:IRGSP-1.0:1:1:43270923:1> rdfs:subClassOf http://www.southgreen.fr/agrold/resource/oryza_sativa/IRGSP-1.0/chromosome:IRGSP-1.0:1:1:43270923:1> .
                    # in taxon = http://purl.obolibrary.org/obo/RO_0002162 ; change taxon	taxon	ObjectProperty	SIO_000253	exact match (SIO_000253:has source) to http://purl.obolibrary.org/obo/RO_0002162
                    # <http://identifiers.org/ensembl/Os12t0534500-00-E1> rdf:type identifiers:ensembl
                    # rdfs:seeAlso panther: <http://purl.uniprot.org/panther/>
                    print(line)
                    ttl_handle.write(line)
    else:
        print("***************** File not found ********************\n")
        print(files)


if len(sys.argv):
    filepath = sys.argv[1]
    ROOT_DIR = os.path.dirname(filepath)
    # print(ROOT_DIR)
    file_input = os.path.basename(filepath)
    # print(file_input)
else:
    ROOT_DIR = '/Volumes/LaCie/AGROLD/data_update_2019/ensembl/'
    file_input = 'zea_mays_xrefs.ttl'

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