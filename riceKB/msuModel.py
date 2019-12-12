from __future__ import print_function
import sys

print(sys.path)
from riceKB.globalVars import *
from riceKB.gffParser import *
import pprint
import re
import os
import pandas as pd
import numpy as np
'''
Created on May, 2017
Updated on Dec, 2019
The msuModel module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling MSU data
It runs with several files downloaded from MSU annotation project http://rice.plantbiology.msu.edu/
all.locus_brief_info.7.0.txt, all.interpro.txt, all.pfam.txt, all.GOSlim_assignment.txt
1 - First define your output path to filename 
path_output = 'all.locus_brief_info.7.ttl'
2 - Pass the arguments to genParser function
ds = geneParser('all.locus_brief_info.7.0.txt',\
                'all.interpro.txt',\
                'all.pfam.txt',\
                'all.GOSlim_assignment.txt')

3 - run the program 
msuModeleRDF(ds, path_output)

@author: larmande
'''
__author__  = "larmande"

#TODO better Error handling


def geneParser(infile, interpro, pfam, go):
    #    pp = pprint.PrettyPrinter(indent=4)
    gene_hash = {}
    tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.\d$')
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    array = pd.read_csv(infile, sep="\t", delimiter=None, dtype='str')
    interpro_array = pd.read_csv(interpro, sep="\t", delimiter=None, header=None,dtype='str')
    pfam_array = pd.read_csv(pfam, sep="\t", delimiter=None, header=None, dtype='str')
    go_array = pd.read_csv(go, sep="\t", delimiter=None, header=None,dtype='str')

    #array['locus_id'].replace('', np.nan, inplace=True)
    #array.dropna(subset=['locus_id'], inplace=True)
    return array, interpro_array, pfam_array, go_array

def msuModeleRDF(msu_ds, output_file):
  # The differentes variable declaration
    line_number = 0
    rdf_writer = open(output_file, "w")
    chromosome_list = list()
    gene_list = list()
    gene,inter,pfam,go = msu_ds
    taxon_id="39947"
    chromosome_size = {'39947': {'size': ['1:1-43270923:1', '2:1-35937250:1', '3:1-36413819:1', '4:1-35502694:1',
                                      '5:1-29958434:1', '6:1-31248787:1', '7:1-29697621:1', '8:1-28443022:1',
                                      '9:1-23012720:1', '10:1-23207287:1', '11:1-29021106:1', '12:1-27531856:1',
                                      'Mt:1-402710:1', 'Pt:1-134481:1'],
                             'genome_assembly': 'IRGSP-1.0'}}
# The first wrinting in the file is the prefix


    print ("*************msu RDF conversion begins***********\n")
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    rdf_writer.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n")
    rdf_writer.write(pr + "\t" + interpro_ns + "<" + interpro_uri + "> .\n")
    rdf_writer.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n")

    # Ajout du prefix pour la realese des donnees
    #rdf_writer.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")

# In here we buil the modele and writer in file with ttl format
    os_japonica_buffer = ''
    #os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + rdfs_ns + "subClassOf" + "\t\t" + sio_ns + "SIO_000253" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "OBI_0100026" + " .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + skos_ns + "prefLabel" + "\t\t" + "\"Oryza sativa Japonica Group\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + rdfs_ns + "label" + "\t\t" + "\"Oryza sativa Japonica Group\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + skos_ns + "altLabel" + "\t\t" + "\"Japanese rice\"" + "@en .\n"
    os_japonica_buffer += ncbi_tax_ns + taxon_id + "\t\t" + dc_ns + "identifier" + "\t\t" + taxon_id + " .\n\n"
    # os_japonica_buffer += ncbi_tax_ns + "39947" + "\t\t" + base_vocab_ns + "taxon" + "\t\t" + ncbi_tax_ns + "39947" + " .\n\n"
    print(os_japonica_buffer)
    rdf_writer.write(os_japonica_buffer)

    for records in gene.as_matrix(columns=None):
        line_number+=1
        # Chromosome triple
        if not records[0] in chromosome_list:
            os_japonica_buffer = ''
            chromosome_list.append(records[0])
            os_japonica_buffer += chromosome_ns + re.sub('Os|Chr', '', records[0]) + "\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "inAssembly" + "\t" +  "\"Os-Nipponbare-Reference-IRGSP-1.0\"" + " ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "inSchemaNumber" + "\t"  + "\"7\"" + " ;\n"
           # os_japonica_buffer += "\t" + base_vocab_ns + "type" + "\t" + res_ns + "Chromosome" + " ;\n\n"
            os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
            rdf_writer.write(os_japonica_buffer)
            print(os_japonica_buffer)

        if len(records) is  10:
            os_japonica_buffer = ''
            if not records[1] in gene_list:
                # print the corresponding gene associated at mRNAs
                #os_japonica_buffer = ''
                gene_list.append(records[1])
                os_japonica_buffer += base_resource_ns + records[1] + "\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + 'IRGSP-1.0' + "\" ;\n"
                os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
                # os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
                os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records[1] + "\" ;\n"
                # os_japonica_buffer += "\t" + rdfs_ns + "subClassOf" + "\t\t" + obo_ns + "SO_0000704" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + " \"" + records[9] + "\" ;\n"
                os_japonica_buffer += "\t" + obo_ns + "obo:RO_0002162" + "\t\t" + ncbi_tax_ns + "39947" + " ;\n"
                os_japonica_buffer += "\t" + base_vocab_ns + "isLocatedOn" + "\t\t" + "" + chromosome_ns + re.sub(
                    'Os|Chr',
                    '',
                    records[
                        0]) + " ;\n"
                os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)
                #rdf_writer.write(os_japonica_buffer)
                #print(os_japonica_buffer)
            # extact on mRNA ids and associated triples
            os_japonica_buffer += base_resource_ns + records[2] + "\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + 'IRGSP-1.0' + "\" ;\n"
            os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
            # os_japonica_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            os_japonica_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records[2] + "\" ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t\t" + base_resource_ns + records[1] + " ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + " \"" + records[9] + "\" ;\n"
            os_japonica_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns  + "39947" + " ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "hasStartPosition" + "\t" + " \"" + str(
            records[3]) + "\"^^xsd:integer ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "hasEndPosition" + "\t" + " \"" + str(
            records[4]) + "\"^^xsd:integer ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "isLocatedOn" + "\t\t" + "" + chromosome_ns + re.sub(
                'Os|Chr',
                '',
                records[
                    0]) + " ;\n"
            os_japonica_buffer += "\t" + base_vocab_ns + "strand" + "\t" + " \"" + records[5] + "\"^^xsd:string ;\n"
            if records[6] == "Y":
                os_japonica_buffer += "\t" + base_vocab_ns + "isTE" + "\t" + "\"true\"^^xsd:boolean ;\n"
            else:
                os_japonica_buffer += "\t" + base_vocab_ns + "isTE" + "\t" + "\"false\"^^xsd:boolean ;\n"

            if records[7] == "Y":
                os_japonica_buffer += "\t" + base_vocab_ns + "isExpressed" + "\t" + "\"true\"^^xsd:boolean ;\n"
            else:
                os_japonica_buffer += "\t" + base_vocab_ns + "isExpressed" + "\t" + "\"false\"^^xsd:boolean ;\n"

            if records[8] == "Y":
                os_japonica_buffer += "\t" + base_vocab_ns + "isRepresentative" + "\t" + "\"true\"^^xsd:boolean ;\n"
            else:
                os_japonica_buffer += "\t" + base_vocab_ns + "isRepresentative" + "\t" + "\"false\"^^xsd:boolean ;\n"
            # extract and print GO ids from go_array
            df_go= go.loc[go[0] == records[2]]
            for index, row in df_go.iterrows():
                if row[2]=="F":
                    os_japonica_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + re.sub(':', '_',\
                                                                                                row[1]) + " ;\n"
                elif row[2]=="P":
                    os_japonica_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + re.sub(':', '_',\
                                                                                                row[1]) + " ;\n"
                elif row[2]=="C":
                    os_japonica_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + re.sub(':', '_', \

                                                                                               row[1]) + " ;\n"
            df_interpro = inter.loc[inter[0] == records[2]]
            lookup_list = list()
            # extract and print interpro ids from interpro_array
            for index, row in df_interpro.iterrows():
                if interpro_pattern.match(str(row[7])) and str(row[7]) not in lookup_list:
                    lookup_list.append(str(row[7]))
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + "<" + up_base_uri + "interpro" + "/" + str(row[7]) + ">" + " ;\n"
            df_pfam = pfam.loc[inter[0] == records[2]]
            # extract and print pfam ids from pfam_array
            lookup_list2 = list()
            for index, row in df_pfam.iterrows(): # use itertuples() instead
                pfam_id = row[1].split('.')[0]
                if pfam_pattern.match(str(pfam_id)) and pfam_id not in lookup_list2:
                    lookup_list2.append(pfam_id)
                    os_japonica_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + "<" + up_base_uri + "pfam" + "/" + str(pfam_id) + ">" + " ;\n"
            # replace the last . to ; triples
            os_japonica_buffer = re.sub(' ;$', ' .\n', os_japonica_buffer)

            rdf_writer.write(os_japonica_buffer)
            print(os_japonica_buffer)
            # to test remove the comment bellow
            # break
    print(line_number)



pp = pprint.PrettyPrinter(indent=4)

#TEST PARAM
#path = '/Users/plarmande/Downloads/IRGSP-1.0_representative/transcripts_mRNA.gff'
path_output = '/Users/plarmande/Downloads/all.locus_brief_info.7.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
#ds = parseGFF3(path)   # The parsing file withe tropGeneParser()
#pp.pprint(ds)    # For to see in teminal the parsing

#ds = os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)
ds = geneParser('/Users/plarmande/Downloads/all.locus_brief_info.7.0.txt',\
                '/Users/plarmande/Downloads/all.interpro.txt',\
                '/Users/plarmande/Downloads/all.pfam.txt',\
                '/Users/plarmande/Downloads/all.GOSlim_assignment.txt')
msuModeleRDF(ds, path_output)
