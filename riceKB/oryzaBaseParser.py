#!/usr/bin/env python
'''
Created on July 17, 2014
The oryzaBaseParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling OryzaBase data

#TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: Larmande
'''
import pprint
from riceKB.globalVars import *
from riceKB.utils import *
import re
import os, sys
import datetime
import json
import pandas as pd
import numpy as np
import rdflib
from rdflib.graph import Graph

'''
OryzaBase Fields
        trait_gene_id, 0
        symbol, 1 list (Note!! entries separated by ',', some entries have '[]', '/', '#')
        alternative_name, 2 list
        name_en, 3 
        allele, 4 list
        chromosome_no, 5
        RAP_id, 6 list
        is_mutant, 7
        arm, 8
        locus, 9
        explanation_en, 10
        recommended_gene_symbol, 11 
        recommended_gene_name, 12 (Note!! some entries have '_')
        protein_name, 13 (Note!! some entries have '_')
        Class name en, 14
        Gene Ontology IDs, 15 list
        Trait Ontology IDs, 16 list
        Gramene ID 17 list (Note!! IDs separated by '/' or/and ',')
'''
rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
gramene_pattern = re.compile(r'^GR\:\d{7}$')
prot_pattern = re.compile(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')
def oryzaBaseParser(input_file):
    array = pd.read_csv(input_file, sep='\t', delimiter=None, dtype='str', skip_blank_lines=True)
    array.fillna('', inplace=True)
    oryGene_ds = dict()
    #fileHandle = open(input_file, 'r')
    #lines  = fileHandle.readlines()
    #lines.pop(0)
#    pp = pprint.PrettyPrinter(indent=4)
    for records in array.to_numpy():
        # current_line = re.sub('\n$', '', current_line)
        #current_line = re.sub('\\r|\\n', '', current_line)
        #records = current_line.split('\t')

        oryGeneID = records[0]
#        pp.pprint(records)
#        print oryGeneID + "\t" + records[6] + "\t" + records[7] + "\t" + records[8]
        #TODO split name column with ',' first row -> name and other -> Alt_names
        #TODO split trait_class with ','
        # TODO split explanation column
        # TODO check chromosome number wrong parsing
        # TODO check if label is fill up with CGSNL Gene Symbol	or	CGSNL Gene Name
        oryGene_ds[oryGeneID] = {
                                  "Symbols": [],
                                  "Alt_names": [],
                                  "Name": records[3],
                                  "Alleles": records[6],
                                  "Chromosome": records[7],#4
                                  "RAP_id": [],
                                  "MSU_id": [],
                                  "Mutant": records[6],
                                  "Arm": records[13],#7
                                  "Locus": records[14],#8
                                  "Explanation": records[8],#9
                                  "Reco_symbol": '',
                                  "Reco_name": '',
                                  "Protein_name" : records[5],#12
                                  "Trait_class": [],#13
                                  "GO_id": [],
                                  "TO_id": [],
                                  "PO_id": [],
                                  "Gramene_id": []
                                  }
        if records[2]:
            symbols = str(records[2]).split(',')
            # symbols = [x.strip() for x in symbols]
            oryGene_ds[oryGeneID]["Symbols"].extend(symbols)
        if records[4]:
            str(records[4]).strip()
            alt_names = str(records[4]).split(',')
            # alt_names = [x.strip() for x in alt_names]
            oryGene_ds[oryGeneID]["Alt_names"].extend(alt_names)
        # if records[5]:
        #     alleles = records[5].split(',')
        # #     # alleles = [x.strip() for x in alleles]
        #     oryGene_ds[oryGeneID]["Alleles"] = records[5]
        if records[9]:
            traits = str(records[9]).split(',')
            # traits = [x.strip() for x in traits]
            oryGene_ds[oryGeneID]["Trait_class"].extend(traits)
        if records[10]:
            rapIds = str(records[10]).split(',|\||\/')
            oryGene_ds[oryGeneID]["RAP_id"].extend(rapIds)
        if records[1] and records[1] != "_":
            oryGene_ds[oryGeneID]["Reco_symbol"] = str(records[1])
        if records[3] and records[3] != "_":
            oryGene_ds[oryGeneID]["Reco_name"] = str(records[3])
        if records[15]:
            go_ids = re.split(r',',str(records[15]))
            oryGene_ds[oryGeneID]["GO_id"].extend(go_ids)
        if records[16]:
            to_ids = re.split(r',', str(records[16]))
            oryGene_ds[oryGeneID]["TO_id"].extend(to_ids)
        if records[17]:
            po_ids = re.split(r',',  str(records[17]))
            oryGene_ds[oryGeneID]["PO_id"].extend(po_ids)
        if records[11]:
            records[11] = re.sub('\\r$', '', str(records[11]))
            ids = re.split(r',|\||\/', str(records[11]))
            # gr_ids = [x.rstrip() for x in gr_ids]
            # oryGene_ds[oryGeneID]["Gramene_id"].extend(gr_ids)
            oryGene_ds[oryGeneID]["MSU_id"].extend(ids)
        if records[12]:
            records[12] = re.sub('\\r$', '', str(records[12]))
            gr_ids = re.split(r'/|,', str(records[12]))
            gr_ids = [x.rstrip() for x in gr_ids]
            oryGene_ds[oryGeneID]["Gramene_id"].extend(gr_ids)
            #oryGene_ds[oryGeneID]["Gramene_id"] = records[10]
    #fileHandle.close()
    return oryGene_ds            
#    pp.pprint(oryGene_ds)
#        print records[0] + "\t" + records[4] 
    
def oryzaBaseRDF(infile, output_file):
#    pp = pprint.PrettyPrinter(indent=4)
    taxon_id = "39947"
    print("*************Parsing OryzaBase gene data ***********\n")
    
    orygene_ds = oryzaBaseParser(infile)
    gene_count = len(orygene_ds)
    
    print("Number of genes: %s\n" % (str(gene_count)))
    print("OryzaBase gene data has been parsed!\n")
    print("*************************************\n\n")
    
    ttl_handle = open(output_file, "w")
    pub_handle = open(pub_dict,'w')
    pub_handle_rapdb = open(pub_dict_RAPDB,'w')
    pub_handle_msu = open(pub_dict_MSU,'w')
    pub_handle_genomeHub = open(pub_dict_HUB,'w')
    PO_handle_genomeHub = open(dict_PO,'w')
    TO_handle_genomeHub = open(dict_TO,'w')
    ttl_buffer = ''
    
    print("************* OryzaBase RDF conversion begins***********\n")


    ttl_handle.write(str(getRDFHeaders()))
#    pp.pprint(orygene_ds)
    for oryid in orygene_ds:
        rap_id = ''
        msu_id = ''
        ttl_buffer = ''
        pub_buffer = ''
        pub_buffer_rapdb = ''
        ttl_buffer += oryzabase_ns + oryid + "\n"
        ttl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
        ttl_buffer += "\t" + dcterms_ns + "identifier" + "\t" + '"%s"' % (oryid) + " ;\n"
        # ttl_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
        for item in orygene_ds[oryid]:
            if item == 'Reco_symbol':
                if orygene_ds[oryid][item]:
                    label =  re.sub('^\s+|^\t+|\s+$', '', orygene_ds[oryid][item])
                    label = re.sub('\"|\'','', label)
                    label = re.sub(r"\\",'', label)
                    ttl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % ( label ) + " ;\n"
                    ttl_buffer += "\t" + skos_ns + "prefLabel" + "\t" + '"%s"' % ( label) + " ;\n"
                    if label != '_':
                        pub_buffer += label  + "\t" + oryzabase_uri + oryid + "\n"
                    pubRAPDB_writer(label,orygene_ds,oryid,pub_handle_rapdb)
                    pubMSU_writer(label, orygene_ds, oryid, pub_handle_msu)
            if item == 'Reco_name':
                if orygene_ds[oryid][item]:
                    description = re.sub('\"|\'','', (orygene_ds[oryid][item]))
                    description = re.sub(r"\\",'', description)
                    ttl_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' %  (description)  + " ;\n"
                    # if description != '_':
                        # pub_buffer += description + "\t" + oryzabase_uri + oryid + "\n"
                    # pubRAPDB_writer(description, orygene_ds, oryid, pub_handle_rapdb)
                    # pubMSU_writer(description, orygene_ds, oryid, pub_handle_msu)
            ## TODO : fix error - sometimes list of name separate by , check how to turn one name and alternate names, variable have "" so remove them
            if item == 'Name':
                if orygene_ds[oryid][item]:
                    altlabel = re.sub('\"|\'','', (orygene_ds[oryid][item]))
                    altlabel = re.sub(r"\\",'', altlabel)
                    ttl_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' %  (altlabel) + " ;\n"
                    # if altlabel != '_':
                    #     pub_buffer += altlabel + "\t" + oryzabase_uri + oryid + "\n"
                    # pubRAPDB_writer(altlabel, orygene_ds, oryid, pub_handle_rapdb)
                    # pubMSU_writer(altlabel, orygene_ds, oryid, pub_handle_msu)
            if item == 'Explanation':
                if orygene_ds[oryid][item]:
                    #bad_chars = '|'.join(['\"', '\'', '\\'])
                    comment = re.sub('\"|\'','', orygene_ds[oryid][item])
                    comment = re.sub(r"\\",'', comment)
                    ttl_buffer += "\t" + rdfs_ns + "comment" + "\t" + '"%s"' %  comment + " ;\n"
            if item == 'Chromosome':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "chromosome" + "\t" + '"%s"' % re.sub('\"|\'', '',(str(orygene_ds[oryid][item])))  + " ;\n"
            # if item == 'Alleles':
            #     if orygene_ds[oryid][item]:
            #         for allele in orygene_ds[oryid][item]:
            #             allele = re.sub('^\s|\s$', '', allele)
            #             ttl_buffer += "\t" + base_vocab_ns + "has_allele" + "\t" + '"%s"' % (allele) + " ;\n"
            if item == 'Mutant':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasMutant" + "\t" + '"%s"' %  re.sub('\"|\'', '',(str(orygene_ds[oryid][item])))  + " ;\n"
            if item == 'Arm':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasChromosomeArm" + "\t" + '"%s"' % re.sub('\"|\'', '',(str(orygene_ds[oryid][item]))) + " ;\n"
            if item == 'Locus':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasLocus" + "\t" + '"%s"' % (orygene_ds[oryid][item]) + " ;\n"
            ## TODO : fix empty fields "" test if value is not null ##
            if item == 'Symbols':
                if orygene_ds[oryid][item]:
                    for symbol in orygene_ds[oryid][item]:
                        if symbol is not None:
                            symbol = re.sub(';', '/',(symbol))
                            symbol = re.sub('^\s+|^\t+|\n+|\s+$', '', symbol)
                            symbol = re.sub(r"\\", '', symbol)
                            if symbol != '_':
                                ttl_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % re.sub('\"|\'', '', symbol )+ " ;\n"
                                pub_buffer += symbol + "\t" + oryzabase_uri + oryid + "\n"
                                pubRAPDB_writer(symbol, orygene_ds, oryid, pub_handle_rapdb)
                                pubMSU_writer(symbol, orygene_ds, oryid, pub_handle_msu)
            if item == 'Alt_names':
                if orygene_ds[oryid][item]:
                    for alt_name in orygene_ds[oryid][item]:
                        alt_name = re.sub('^\s+|^\t+|\n+|\s+$', '', alt_name)
                        alt_name = re.sub('\"+', '', alt_name)
                        alt_name = re.sub(r"\\",'', alt_name)
                        if alt_name != '_':
                            ttl_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (alt_name) + " ;\n"
                            # pub_buffer += alt_name + "\t" + oryzabase_uri + oryid + "\n"
                            # pubRAPDB_writer(alt_name, orygene_ds, oryid, pub_handle_rapdb)
                            # pubMSU_writer(alt_name, orygene_ds, oryid, pub_handle_msu)
            if item == 'RAP_id':
                if orygene_ds[oryid][item]:
                    for rap_id in orygene_ds[oryid][item]:
                        if re.match(rap_pattern, rap_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + re.sub('\s+', '', rap_id)  + " ;\n"
                            ttl_buffer += "\t" + owl_ns + "sameAs" + "\t" + res_ns + re.sub('\s+', '', rap_id)  + " ;\n"
                            pub_buffer += rap_id + "\t" + rapdb_gene_uri + re.sub('\s+', '', rap_id) +"\n"
                            pub_handle_genomeHub.write(re.sub('\s+', '', rap_id) + "\t" + ''.join(orygene_ds[oryid]['Reco_symbol']) + "\t" + ''.join(orygene_ds[oryid]['Symbols']) + "\n")
            if item == 'MSU_id':
                if orygene_ds[oryid][item]:
                    for msu_id in orygene_ds[oryid][item]:
                        if re.match(tigr_pattern, msu_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + res_ns + re.sub('\s+', '',
                                                                                               msu_id) + " ;\n"
                            ttl_buffer += "\t" + owl_ns + "sameAs" + "\t" + res_ns + re.sub('\s+', '', msu_id) + " ;\n"
                            pub_buffer += msu_id + "\t" + sniplay_gene_uri + re.sub('\.\d\s*$','',msu_id) + "\n"
            if item == 'Gramene_id':
                if orygene_ds[oryid][item]:
                    for gr_id in orygene_ds[oryid][item]:
                        if re.match(gramene_pattern,gr_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + gr_g_ns + re.sub('\s+', '',gr_id) + " ;\n"

            if item == 'Protein_name':
                if orygene_ds[oryid][item]:
                    if re.match(prot_pattern, orygene_ds[oryid][item]):
                        ttl_buffer += "\t" + sio_ns + "SIO_010078" + "\t" + '"%s"' % re.sub('\"|\'', '',(str(orygene_ds[oryid][item]))) + " ;\n" # encodes
            if item == 'Trait_class':
                if orygene_ds[oryid][item]:
                    for traits in orygene_ds[oryid][item]:
                         for trait in traits.split('-'):
                            trait = re.sub('^\s+|\s+$', '', trait)
                            ttl_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + '"%s"' %  re.sub('\"|\'', '',(trait))+ " ;\n"
            if item == 'TO_id':
                if orygene_ds[oryid][item]:
                    for to_term in orygene_ds[oryid][item]:
                        for to_id in to_term.split('-'):
                            to_id = re.sub(':', '_', to_id)
                            to_id = re.sub('^\s|\s$', '', to_id)
                            if re.match(r"^TO_[0-9]{7}$", to_id):
                                ttl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + obo_ns + to_id + " ;\n"
                        if rap_id is not None and rap_id != "":
                            TO_handle_genomeHub.write(rap_id + "\t" + re.sub('^\s|\s$', '', to_term) + "\n")
                        # if re.match(r"^TO_[0-9]{7} or TO_[0-9]{7}", to_term):
                        #     TO_oryza = to_term.split(' or ')
                        #     for TO in TO_oryza:
                        #         print(TO)
                        #     ttl_buffer += "\t" + base_vocab_ns + "has_trait" + "\t" + obo_ns + TO + " ;\n"
                        # else:
                        #     ttl_buffer += "\t" + base_vocab_ns + "has_trait" + "\t" + '"' + to_term + '"' + ";\n"

            if item == 'GO_id':
                if orygene_ds[oryid][item]:
                    for go_term in orygene_ds[oryid][item]:
                        for go_id in go_term.split('-'):
                            go_id = re.sub(':', '_', go_id)
                            go_id = re.sub('\s', '', go_id)
                            if re.match(r"^GO_[0-9]{7}$", go_id):
                                ttl_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + go_id + " ;\n"
            if item == 'PO_id':
                if orygene_ds[oryid][item]:
                    for po_term in orygene_ds[oryid][item]:
                        for po_id in po_term.split('-'):
                            po_id = re.sub(':', '_', po_id)
                            po_id = re.sub('\s', '', po_id)
                            if re.match(r"^TO_[0-9]{7} or TO_[0-9]{7}", po_id):
                                ttl_buffer += "\t" + obo_ns + "RO_0002206" + "\t" + obo_ns + po_id + " ;\n" # expressed in
                        if rap_id is not None and rap_id != "":
                            PO_handle_genomeHub.write(rap_id + "\t" + re.sub('^\s|\s$', '', po_term) + "\n")
            
        ttl_buffer = re.sub(' ;$', ' .\n', ttl_buffer)
        pub_handle.write(pub_buffer)
        #pub_handle_rapdb.write(pub_buffer_rapdb)
        RDF_validation(ttl_buffer,ttl_handle,oryid)
        #ttl_handle.write(ttl_buffer)

    ttl_handle.close()
    print("************* OryzaBase RDF completed ************!\n\n")

    print("************* OryzaBase RDF Testing Graph ************!\n\n")
def pubRAPDB_writer(label,orygene_ds,oryid,pub_handle_rapdb):
    if orygene_ds[oryid]['RAP_id']:
        for rap_id in orygene_ds[oryid]['RAP_id']:
            rap_id = re.sub('^\s+|^\t+|\s+$', '', rap_id)
            if rap_id.split('/'):
                rap_dic = rap_id.split('/')
                for rap_id in rap_dic:
                    if re.match(rap_pattern, rap_id):
                        if label != '_':
                            pub_handle_rapdb.write(label + "\t" + rapdb_gene_uri + rap_id + "\n")
                        #print(rap_id)
            else:
                if re.match(rap_pattern, rap_id):
                    if label != '_':
                        pub_handle_rapdb.write(label + "\t" + rapdb_gene_uri + rap_id + "\n")

def pubMSU_writer(label, orygene_ds, oryid, pub_handle_msu):
    if orygene_ds[oryid]['MSU_id']:
        for msu_id in orygene_ds[oryid]['MSU_id']:
            msu_id = re.sub('^\s+|^\t+|\s+$', '', msu_id)
            msu_id = re.sub('\.\d\s*$','',msu_id)
            if re.match(tigr_pattern, msu_id):
                if label != '_':
                    pub_handle_msu.write(label + "\t" + sniplay_gene_uri + msu_id + "\n")

# def dicTO_writer(label, orygene_ds, oryid, pub_handle_msu):

def RDF_validation(ttl_buffer,ttl_handle,oryid):

    try:
        temp_file = '/Users/pierre/Downloads/tmp/temp_graph.ttl'
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
        temp = '/Users/pierre/Downloads/tmp/temp_graph'+ oryid +'.ttl'
        handle = open(temp, "w")
        handle.write(str(getRDFHeaders()))
        handle.write(ttl_buffer)
        pass

oryzabase_file = '/Users/pierre/workspace2015/datasets/OryzabaseGeneListEn_20221122010058.txt'


oryzaBase_output = '/Users/pierre/workspace2015/datasets/OryzabaseGeneListEn_20221122010058.ttl'
pub_dict = '/Users/pierre/workspace2015/datasets/pub_dictionnary.txt'
pub_dict_RAPDB = '/Users/pierre/workspace2015/datasets/pub_dictionnary_rapdb.txt'
pub_dict_MSU = '/Users/pierre/workspace2015/datasets/pub_dictionnary_msu.txt'
pub_dict_HUB = '/Users/pierre/workspace2015/datasets/pub_genome_hub.txt'
dict_PO = '/Users/pierre/workspace2015/datasets/dict_PO_hub.txt'
dict_TO = '/Users/pierre/workspace2015/datasets/dict_TO_hub.txt'
oryzaBaseRDF(oryzabase_file, oryzaBase_output)

#
# g = Graph()
# g.parse(oryzaBase_output, format="turtle")
#
# print len(g)