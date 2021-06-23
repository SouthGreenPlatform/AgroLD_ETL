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
Created on Jan, 2021
This module is created as part of the AgroLD Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Pearl Millet genome annotation \
datasets to normalize with AgroLD

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

def milletParser(files, type):
    """
       parse the ensembl rdf file
       :param files: path of the file
       write  objects corresponding to prefix and triples
    """
    if (type == "annotation"):
        print('OK annotation')

        file1 = "annotation_modified.tsv"
        file_path = os.path.join(data_dir, file1)
        annotation_hash = {}

        with open(file_path, "r") as filehandle:
            for line in filehandle:
                line = line.rstrip('\n')
                line = re.sub('\"', '', line)

                records = line.split('\t')
                geneid = records.pop(0)
                db = records.pop(0)
                if geneid in annotation_hash:
                    annotation_hash[geneid].append({'type': db, 'value': records})
                else:
                    annotation_hash[geneid]=[]
                    annotation_hash[geneid].append({'type':db, 'value': records})
        return annotation_hash
    else:

        file1 = "liste_gene_autre.txt"
        path = os.path.join(data_dir, file1)
    #output_writer = open(output_file, "w")
        print("*************Parsing Pearl Millet gene list data ***********\n")
        array = pd.read_csv(path, sep='\t', delimiter=None, dtype='str',\
                                              skip_blank_lines=True)
        array.fillna('', inplace=True)
        print(array.head())
        row_count = array.shape[0]
        print("Number of genes: %s\n" % (str(row_count)))
        print(" gene list data has been parsed!\n")
        print("*************************************\n\n")
        return array

def gene2RDF(gene,output):
    print("************* Gene RDF conversion begins***********\n")
    ttl_handle = open(path_output, "w")
    ttl_handle.write(str(getRDFHeaders()))
    for records in gene.to_numpy():
        geneid = records[0]
        (strand, position) = getStrandValue(records[3])
        taxon_id = "4543"
        rdf_buffer = ''
        rdf_buffer += "<" + base_resource_uri +  geneid + ">\n"
        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + obo_ns + "SO_0000704" + " ;\n"
        rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + geneid + "\" ;\n"
        rdf_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        rdf_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + geneid + "\" ;\n"
        chromosome_number = re.sub('chr','', str(records[4]))
        rdf_buffer += "\t" + faldo_ns + "location" + "\t" + "<http://www.southgreen.fr/agrold/resource/chromosome/" + taxon_id + "/" + \
                    chromosome_number + ':' + str(records[1]) + '-' + str(records[2]) + ":" + records[3] + ">  .\n\n"

        rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
        rdf_buffer += getFaldoRegion(taxon_id,chromosome_number,records[1],records[2],records[3])
        # Region
        # rdf_buffer += chromosome_ns + taxon_id + "/" + \
        #                       chromosome_number + ':' + str(records[1]) + '-' + str(records[2]) + ":" + strand + "  \n"
        # rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + chromosome_ns + taxon_id + "/" + \
        #                       chromosome_number + ':' + str(records[1]) + '-' + str(records[2]) + ":" + strand + "\";\n"
        # rdf_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
        # rdf_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + taxon_id + "/" + + \
        #                       chromosome_number + ":" + str(records[1]) + ":" + strand + "  ;\n"
        # rdf_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + taxon_id + "/" + + \
        #                       chromosome_number + ":" +  str(records[2]) + ":" + strand + "  .\n\n"
        #
        # # Position 1
        # rdf_buffer += chromosome_ns + "CEGSBv1.1:" + chromosome_number + ":" + str(records[1]) + ":" + strand
        # rdf_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
        # rdf_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
        # rdf_buffer += "  ;\n"
        # rdf_buffer += "\t" + faldo_ns + "position" + "\t" + str(records[1]) + " ;\n"
        # rdf_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "CEGSBv1.1:" + chromosome_number \
        #               +  " .\n\n"
        #
        # # Position 2
        # rdf_buffer += chromosome_ns + "CEGSBv1.1:" + chromosome_number + ":" + str(records[2]) + ":" + strand
        # rdf_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
        # rdf_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
        # rdf_buffer += "  ;\n"
        # rdf_buffer += "\t" + faldo_ns + "position" + "\t" + str(records[2]) + " ;\n"
        # rdf_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + "CEGSBv1.1:" \
        #               + chromosome_number +  " .\n\n"

        #print(rdf_buffer)
        # rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
        #RDF_validation(rdf_buffer, ttl_handle, geneid)
        ttl_handle.write(rdf_buffer)
    ttl_handle.close()
def annotation2RDF(annotation, output):
    print("************* Annotation RDF conversion begins***********\n")
    count = 0
    ttl_handle = open(path_output, "a")
    prot_handle = open(uniprotid_list, 'w')
    prot_buffer = ''
    #ttl_handle.write(str(getRDFHeaders()))

    for geneid in annotation:
        annotation_list = annotation[geneid]
        count += 1
        protlist = set()
        taxon_id = "4543"
        rdf_buffer = ''
        rdf_buffer += "<" + base_resource_uri + "transcript/" + geneid + ">\n"
        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
        rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + geneid + "\" ;\n"
        rdf_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
        rdf_buffer += "\t" + dcterms_ns + "identifier" + "\t" + " \"" + geneid + "\" ;\n"
        # rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t\t" + ensembl_transcript_ns + records[0] + ";\n"
        rdf_buffer += "\t" + sio_ns + "SIO_010081" + "\t\t" + res_ns + geneid + " ;\n"
        # rdf_buffer += "\t" + dc_ns + "description" + "\t" + "\"%s" % (
        #     re.sub('\"|\'', '', str(records[2]))) + "\" ;\n"
        #print(len(annotation_list))
        for i in range(len(annotation_list)):
            annotation_dict = annotation[geneid][i]

            if (annotation_dict['type']== 'GO'):
                for go_tuple in annotation_dict['value']:
                    if go_tuple != '':
                        go_list = go_tuple.split(';')
                        if go_list[0] != '':
                            go_id = re.sub(':', '_', go_list[0])
                            go_id = re.sub('\s', '', go_id)
                            if re.match(r"^GO_[0-9]{7}$", go_id):
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + go_id + " ;\n"
            elif (annotation_dict['type']== 'InterPro'):
                for ipr_tuple in annotation_dict['value']:
                    if ipr_tuple != '':
                        ipr_list = ipr_tuple.split(';')
                        if ipr_list[0] != '':
                            ipr_id = re.sub('\s', '', ipr_list[0])
                            if re.match(interpro_pattern, ipr_id):
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + interpro_ns + ipr_id + " ;\n"
            elif (annotation_dict['type']== 'KEGG'):
                for kegg_tuple in annotation_dict['value']:
                    if kegg_tuple != '':
                        kegg_tuple = re.sub('\r', '', kegg_tuple)
                        rdf_buffer += "\t" + dcterms_ns + "description" + "\t" '"%s"' % kegg_tuple + " ;\n"
            elif (annotation_dict['type']== 'Swiss-Prot'):
                for prot_tuple in annotation_dict['value']:
                    if prot_tuple != '':
                        uniprotid = prot_tuple.split(',')[0]
                        if re.match(prot_pattern, uniprotid):
                            rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + uniprot_ns + uniprotid + " ;\n"
                            rdf_buffer += "\t" + base_vocab_ns + "SEQUENCE_MATCH" + "\t" + uniprot_ns + uniprotid + " ;\n"
                            protlist.add(uniprotid)
            elif (annotation_dict['type'] == 'TrEMBL'):
                for trembl_tuple in annotation_dict['value']:
                    if trembl_tuple != '':
                        tremblid = trembl_tuple.split(',')[0]
                        if re.match(prot_pattern, tremblid):
                            rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + uniprot_ns + tremblid + " ;\n"
                            rdf_buffer += "\t" + base_vocab_ns + "SEQUENCE_MATCH" + "\t" + uniprot_ns + tremblid + " ;\n"
                            protlist.add(tremblid)

        rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
        for prot in protlist:
            #print(prot)
            prot_handle.write(prot+"\n")
        #RDF_validation(rdf_buffer, ttl_handle, geneid)
        ttl_handle.write(rdf_buffer)
    ttl_handle.close()
    prot_handle.close()
#TEST PARAM

#path = '/Users/pierre/workspace2015/datasets/Pearl_Millet/liste_gene_autre.txt'
rdf_file = 'annotation.ttl' # The output
#pub_dict = '/Users/plarmande/workspace2015/datasets/uniprot-dictionnary-millet.txt'
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output]
data_dir = sys.argv.pop() # path to the dataset
ROOT_DIR = sys.argv.pop() # path to the root folder
uniprotid_list = sys.argv.pop() # path
path_output = os.path.join(ROOT_DIR, rdf_file)
print("%s .... %s ...%s .... %s" % (data_dir,ROOT_DIR,uniprotid_list,path_output))

gene = milletParser(data_dir,"gene")
gene2RDF(gene,path_output)
annotation = milletParser(data_dir,"annotation")
#print(annotation)
annotation2RDF(annotation,path_output)