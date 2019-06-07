#!/usr/bin/env python
'''
Created on July 17, 2014
The oryzaBaseParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling OryzaBase data

#TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: venkatesan
'''
import pprint
from globalVars import *
import re
import os, sys
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
    oryGene_ds = dict()
    fileHandle = open(input_file, "rU")
    lines  = fileHandle.readlines()
    lines.pop(0)
#    pp = pprint.PrettyPrinter(indent=4)
    for current_line in lines:
        # current_line = re.sub('\n$', '', current_line)
        current_line = re.sub('\\r|\\n', '', current_line)
        records = current_line.split('\t')

        oryGeneID = records.pop(0)
#        pp.pprint(records)
#        print oryGeneID + "\t" + records[6] + "\t" + records[7] + "\t" + records[8]
        #TODO split name column with ',' first row -> name and other -> Alt_names
        #TODO split trait_class with ','
        # TODO split explanation column
        # TODO check chromosome number wrong parsing
        oryGene_ds[oryGeneID] = {
                                  "Symbols": [],
                                  "Alt_names": [],
                                  "Name": records[2],
                                  "Alleles": records[5],
                                  "Chromosome": records[6],#4
                                  "RAP_id": [],
                                  "Mutant": records[5],
                                  "Arm": records[11],#7
                                  "Locus": records[12],#8
                                  "Explanation": records[7],#9
                                  "Reco_symbol": '',
                                  "Reco_name": '',
                                  "Protein_name" : records[4],#12
                                  "Trait_class": [],#13
                                  "GO_id": [],
                                  "TO_id": [],
                                  "PO_id": [],
                                  "Gramene_id": records[10]
                                  }
        if records[1]:
            symbols = records[1].split(',')
            # symbols = [x.strip() for x in symbols]
            oryGene_ds[oryGeneID]["Symbols"].extend(symbols)
        if records[3]:
            records[3].strip()
            alt_names = records[3].split(',')
            # alt_names = [x.strip() for x in alt_names]
            oryGene_ds[oryGeneID]["Alt_names"].extend(alt_names)
        # if records[5]:
        #     alleles = records[5].split(',')
        # #     # alleles = [x.strip() for x in alleles]
        #     oryGene_ds[oryGeneID]["Alleles"] = records[5]
        if records[8]:
            traits = records[8].split(',')
            # traits = [x.strip() for x in traits]
            oryGene_ds[oryGeneID]["Trait_class"].extend(traits)
        if records[9]:
            rapIds = records[9].split(',|\||\/')
            oryGene_ds[oryGeneID]["RAP_id"].extend(rapIds)
        if records[0] and records[0] != "_":
            oryGene_ds[oryGeneID]["Reco_symbol"] = records[0]
        if records[2] and records[2] != "_":
            oryGene_ds[oryGeneID]["Reco_name"] = records[2]
        if records[13]:
            go_ids = re.split(r',',records[13])
            oryGene_ds[oryGeneID]["GO_id"].extend(go_ids)
        if records[14]:
            to_ids = re.split(r',', records[14])
            oryGene_ds[oryGeneID]["TO_id"].extend(to_ids)
        if records[15]:
            po_ids = re.split(r',',  records[15])
            oryGene_ds[oryGeneID]["PO_id"].extend(po_ids)
        if records[10]:
            records[10] = re.sub('\\r$', '', records[10])
            # gr_ids = re.split(r'/|,', records[10])
            # gr_ids = [x.rstrip() for x in gr_ids]
            # oryGene_ds[oryGeneID]["Gramene_id"].extend(gr_ids)
            oryGene_ds[oryGeneID]["Gramene_id"] = records[10]
    fileHandle.close() 
    return oryGene_ds            
#    pp.pprint(oryGene_ds)
#        print records[0] + "\t" + records[4] 
    
def oryzaBaseRDF(infile, output_file):
#    pp = pprint.PrettyPrinter(indent=4)

    print "*************Parsing OryzaBase gene data ***********\n" 
    
    orygene_ds = oryzaBaseParser(infile)
    gene_count = len(orygene_ds)
    
    print "Number of genes: %s\n" % (str(gene_count))
    print "OryzaBase gene data has been parsed!\n"
    print "*************************************\n\n"
    
    ttl_handle = open(output_file, "w")
    ttl_buffer = ''
    
    print "************* OryzaBase RDF conversion begins***********\n" 
    
    ttl_handle.write(base + "\t" + "<" + base_uri + "> .\n")    
    ttl_handle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    ttl_handle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    ttl_handle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    ttl_handle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    ttl_handle.write(pr + "\t" + gr_g_ns + "<" + gramene_gene + "> .\n")
    ttl_handle.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    ttl_handle.write(pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n\n") #tigr_g_ns tigr_g_uri
#    pp.pprint(orygene_ds)
    for oryid in orygene_ds:
        ttl_buffer = ''
        
        ttl_buffer += oryzabase_ns + oryid + "\n"
        ttl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
        # ttl_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
        for item in orygene_ds[oryid]:
            if item == 'Reco_symbol':
                if orygene_ds[oryid][item]:
                    label =  re.sub('^\s+|^\t+|\s+$', '', orygene_ds[oryid][item])
                    ttl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % re.sub('\"|\'', '', label ) + " ;\n"
            if item == 'Reco_name':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "description" + "\t" + '"%s"' %  re.sub('\"|\'', '',(orygene_ds[oryid][item]))  + " ;\n"
            ## TODO : fix error - sometimes list of name separate by , check how to turn one name and alternate names, variable have "" so remove them
            if item == 'Name':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "name" + "\t" + '"%s"' %  re.sub('\"|\'', '',(orygene_ds[oryid][item])) + " ;\n"
            if item == 'Explanation':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "explanation" + "\t" + '"%s"' %  re.sub('\"|\'', '',(orygene_ds[oryid][item])) + " ;\n"
            if item == 'Chromosome':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "chromosome" + "\t" + '"%s"' % re.sub('\"|\'', '',(orygene_ds[oryid][item]))  + " ;\n"
            # if item == 'Alleles':
            #     if orygene_ds[oryid][item]:
            #         for allele in orygene_ds[oryid][item]:
            #             allele = re.sub('^\s|\s$', '', allele)
            #             ttl_buffer += "\t" + base_vocab_ns + "has_allele" + "\t" + '"%s"' % (allele) + " ;\n"
            if item == 'Mutant':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "has_mutant" + "\t" + '"%s"' %  re.sub('\"|\'', '',(orygene_ds[oryid][item]))  + " ;\n"
            if item == 'Arm':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "has_chromosome_arm" + "\t" + '"%s"' % re.sub('\"|\'', '',(orygene_ds[oryid][item])) + " ;\n"
            if item == 'Locus':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "has_locus" + "\t" + '"%s"' % (orygene_ds[oryid][item]) + " ;\n"                    
            ## TODO : fix empty fields "" test if value is not null ##
            if item == 'Symbols':
                if orygene_ds[oryid][item]:
                    for symbol in orygene_ds[oryid][item]:
                        if symbol is not None:
                            symbol = re.sub(';', '/',(symbol))
                            symbol = re.sub('^\s+|^\t+|\s+|\s+$', '', symbol)
                            ttl_buffer += "\t" + base_vocab_ns + "has_synonym" + "\t" + '"%s"' % re.sub('\"|\'', '', symbol )+ " ;\n"
            if item == 'Alt_names':
                if orygene_ds[oryid][item]:
                    for alt_name in orygene_ds[oryid][item]:
                        alt_name = re.sub('^\s+|^\t+|\s+$', '', alt_name)
                        alt_name = re.sub('\"+', '', alt_name)
                        # alt_name = re.sub("\\", '', alt_name)
                        ttl_buffer += "\t" + base_vocab_ns + "has_alternative_name" + "\t" + '"%s"' % re.sub('\"|\'', '',(alt_name)) + " ;\n"
            if item == 'RAP_id':
                if orygene_ds[oryid][item]:
                    for rap_id in orygene_ds[oryid][item]:
                        if re.match(rap_pattern, rap_id):
                            ttl_buffer += "\t" + base_vocab_ns + "has_rap_identifier" + "\t" + ensembl_ns + re.sub('\s+', '', rap_id)  + " ;\n"
            if item == 'Gramene_id':
                if orygene_ds[oryid][item]:
                    for gr_id in orygene_ds[oryid][item]:
                        if re.match(gramene_pattern,gr_id):
                            ttl_buffer += "\t" + base_vocab_ns + "has_gramene_identifier" + "\t" + gr_g_ns + re.sub('\s+', '',gr_id) + " ;\n"
            if item == 'Protein_name':
                if orygene_ds[oryid][item]:
                    if re.match(prot_pattern, orygene_ds[oryid][item]):
                        ttl_buffer += "\t" + base_vocab_ns + "has_protein_name" + "\t" + '"%s"' % re.sub('\"|\'',(orygene_ds[oryid][item])) + " ;\n"
            if item == 'Trait_class':
                if orygene_ds[oryid][item]:
                    for traits in orygene_ds[oryid][item]:
                         for trait in traits.split('-'):
                            trait = re.sub('^\s+|\s+$', '', trait)
                            ttl_buffer += "\t" + base_vocab_ns + "has_trait_class" + "\t" + '"%s"' %  re.sub('\"|\'', '',(trait))+ " ;\n"
            if item == 'TO_id':
                if orygene_ds[oryid][item]:
                    for to_term in orygene_ds[oryid][item]:
                        for to_id in to_term.split('-'):
                            to_id = re.sub(':', '_', to_id)
                            to_id = re.sub('^\s|\s$', '', to_id)
                            if re.match(r"^TO_[0-9]{7}$", to_id):
                                ttl_buffer += "\t" + base_vocab_ns + "has_trait" + "\t" + obo_ns + to_id + " ;\n"

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
                                ttl_buffer += "\t" + base_vocab_ns + "go_term" + "\t" + obo_ns + go_id + " ;\n"
            if item == 'PO_id':
                if orygene_ds[oryid][item]:
                    for po_term in orygene_ds[oryid][item]:
                        for po_id in po_term.split('-'):
                            po_id = re.sub(':', '_', po_id)
                            po_id = re.sub('\s', '', po_id)
                            if re.match(r"^TO_[0-9]{7} or TO_[0-9]{7}", po_id):
                                ttl_buffer += "\t" + base_vocab_ns + "expressed_in" + "\t" + obo_ns + po_id + " ;\n"
            
        ttl_buffer = re.sub(' ;$', ' .\n', ttl_buffer)
        RDF_validation(ttl_buffer,ttl_handle,oryid)

    ttl_handle.close()
    print "************* OryzaBase RDF completed ************!\n\n"

    print "************* OryzaBase RDF Testing Graph ************!\n\n"

def RDF_validation(ttl_buffer,ttl_handle,oryid):

    try:
        temp_file = '/Users/plarmande/Downloads/temp_graph.ttl'
        try_handle = open(temp_file, "w")
        try_handle.write(base + "\t" + "<" + base_uri + "> .\n")
        try_handle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
        try_handle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
        try_handle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
        try_handle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
        try_handle.write(pr + "\t" + gr_g_ns + "<" + gramene_gene + "> .\n")
        try_handle.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
        try_handle.write(pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n\n")  # tigr_g_ns tigr_g_uri
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
        temp = '/Users/plarmande/Downloads/temp_graph'+ oryid +'.ttl'
        handle = open(temp, "w")
        handle.write(base + "\t" + "<" + base_uri + "> .\n")
        handle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
        handle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
        handle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
        handle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
        handle.write(pr + "\t" + gr_g_ns + "<" + gramene_gene + "> .\n")
        handle.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
        handle.write(pr + "\t" + oryzabase_ns + "<" + oryzabase_uri + "> .\n\n")  # tigr_g_ns tigr_g_uri
        handle.write(ttl_buffer)
        pass

oryzabase_file = '/Users/plarmande/Downloads/OryzabaseGeneListEn_20190528010057.txt'


oryzaBase_output = '/Users/plarmande/Downloads/OryzabaseGeneListEn_20190528010057.ttl'

oryzaBaseRDF(oryzabase_file, oryzaBase_output)

#
# g = Graph()
# g.parse(oryzaBase_output, format="turtle")
#
# print len(g)