import sys
import argparse
from riceKB.globalVars import *
from riceKB.utils import *
import pprint
import re
import os
import pandas as pd
import numpy as np

'''
Created on Dec, 2019
The interpro module is created as part of the AgroLD project.

This module contains Parsers, RDF converters and generic functions for handling interpro data
It runs with ParentChildTreeFile.txt file downloaded from EBI interpro website

1 - set up the file input arguments "path" for the input file and "path_output" for the output file

2 - run the program 
fileParser(path,path_output)

@author: larmande
'''
__author__  = "larmande"

def fileParser(infile,list, mapping, outfile):
    #array = pd.read_csv(infile, sep="\t", delimiter=None , dtype='str')
    # entry.list has a header row (ENTRY_AC/ENTRY_TYPE/ENTRY_NAME); header=0
    # drops it instead of ingesting it as a bogus entry.
    array1 = pd.read_csv(list, sep="\t", delimiter=None, dtype='str' , header=0, skip_blank_lines=True)
    # interpro2go starts with a handful of "!"-prefixed metadata lines (date,
    # citation, contact...) with no " > " in them; comment='!' drops those so
    # pandas doesn't infer a 1-column format from them and then choke on the
    # first real, 2-column row.
    goarray = pd.read_csv(mapping, sep=" > ", delimiter=None, dtype='str' , header=None, engine='python', skip_blank_lines=True, comment='!')
    godic = {}
    for line in goarray.to_numpy():
        iprtuple = line[0]
        iprid = iprtuple.split(' ')[0]
        iprid = re.sub('InterPro:', '', iprid)
        item = line[1]
        goitem = item.split('; ')[1]
        if iprid in godic:
            godic[iprid].append(goitem)
        else:
            godic[iprid] = []
            godic[iprid].append(goitem)

    rdf_writer = open(outfile, "w")
    rdf_writer.write(str(getRDFHeaders()))
    nb_line=0
    for records in array1.to_numpy():
        iprid = records[0]
        nb_line +=1
        buffer = ''
        buffer += "<" + interpro_ns + iprid + ">\n"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Interpro_domain" + " ;\n"
        buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + records[1] + " ;\n"
        buffer += "\t" + rdfs_ns + "label" + "\t" + "\"" + records[2] + "\" ;\n"
        buffer += "\t" + dc_ns + "identifier" + "\t" + "\"" + iprid + "\" ;\n"
        if iprid in godic:
            for goid in godic[iprid]:
                goid = re.sub(':','_',goid)
                if goid != '':
                    buffer += "\t" + base_vocab_ns + "mappingWith" + "\t" + obo_ns + goid + " ;\n"
        buffer = re.sub(' ;$', ' .\n', buffer)
        rdf_writer.write(buffer)
        print(buffer)
    print("nb de lignes traitees:" + str(nb_line))

    with open(infile, 'r', encoding='utf-8') as file:
        root_item = ""
        current_item = ""
        second_item  = ""
        third_item = ""
        fourth_item = ""
        #rdf_writer = open(outfile, "w")
        nb_line = 0
        #rdf_writer.write(str(getRDFHeaders()))
        for line in file:
            buffer = ''
            if not line.startswith('--'):
                nb_line += 1
                key,value = line.rstrip('\n').rstrip(':').split('::', 1)
                root_item = key
            elif not line.startswith('----'):
                nb_line += 1
                key, value = line.rstrip('\n').rstrip(':').split('::', 1)
                key = re.sub('--','',key)
                second_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + root_item + " ;\n"
            elif not line.startswith('------'):
                nb_line += 1
                key, value = line.rstrip('\n').rstrip(':').split('::', 1)
                key = re.sub('--', '', key)
                third_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + second_item + " ;\n"
            elif not line.startswith('--------'):
                nb_line += 1
                key, value = line.rstrip('\n').rstrip(':').split('::', 1)
                key = re.sub('--', '', key)
                fourth_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + third_item + " ;\n"
            else:
                nb_line += 1
                key, value = line.rstrip('\n').rstrip(':').split('::', 1)
                key = re.sub('--', '', key)
                parent_item = current_item
                current_item = key
                buffer += "<" + interpro_ns + key + ">\n"
                buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + interpro_ns + fourth_item + " ;\n"

            buffer = re.sub(' ;$', ' .\n', buffer)
            rdf_writer.write(buffer)
            print(buffer)
        print("nb de lignes traitees:" + str(nb_line))



def main():
    parser = argparse.ArgumentParser(
        description="Convert InterPro's ParentChildTreeFile, entry list and "
                    "GO mapping to RDF.")
    parser.add_argument("tree_file", help="ParentChildTreeFile.txt (entry hierarchy)")
    parser.add_argument("entry_list", help="entry.list (id, type, label)")
    parser.add_argument("go_mapping", help="interpro2go (id -> GO term mapping)")
    parser.add_argument("output", help="output Turtle file")
    args = parser.parse_args()
    fileParser(args.tree_file, args.entry_list, args.go_mapping, args.output)


if __name__ == "__main__":
    main()