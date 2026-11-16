#!/usr/bin/env python
'''
Created on July 17, 2014
The oryzaBaseParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling OryzaBase data

Input: OryzaBase's "Gene list" export (tab-separated, header row), e.g.
OryzabaseGeneListEn_YYYYMMDDHHMMSS.txt, columns (0-indexed, verified against
a real 2026-10-22 export -- the field list below is accurate; don't trust an
older comment here over this one if the two ever disagree again):
    0  Trait Gene Id
    1  CGSNL Gene Symbol
    2  Gene symbol synonym(s)
    3  CGSNL Gene Name
    4  Gene name synonym(s)
    5  Protein Name
    6  Allele
    7  Chromosome No.
    8  Explanation
    9  Trait Class
    10 RAP ID
    11 MSU ID
    12 Gramene ID
    13 Arm
    14 Locate(cM)
    15 Gene Ontology
    16 Trait Ontology
    17 Plant Ontology

@author: Larmande
'''
import argparse
import os
import re
from riceKB.globalVars import *
from riceKB.utils import *
import pandas as pd

rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
gramene_pattern = re.compile(r'^GR\:\d{7}$')
prot_pattern = re.compile(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')


def oryzaBaseParser(input_file):
    array = pd.read_csv(input_file, sep='\t', dtype='str', skip_blank_lines=True)
    array.fillna('', inplace=True)
    oryGene_ds = dict()
    for records in array.to_numpy():
        oryGeneID = records[0]
        oryGene_ds[oryGeneID] = {
                                  "Symbols": [],
                                  "Alt_names": [],
                                  "Name": records[3],
                                  "Alleles": records[6],
                                  "Chromosome": records[7],
                                  "RAP_id": [],
                                  "MSU_id": [],
                                  "Mutant": records[6],
                                  "Arm": records[13],
                                  "Locus": records[14],
                                  "Explanation": records[8],
                                  "Reco_symbol": '',
                                  "Reco_name": '',
                                  "Protein_name" : records[5],
                                  "Trait_class": [],
                                  "GO_id": [],
                                  "TO_id": [],
                                  "PO_id": [],
                                  "Gramene_id": []
                                  }
        if records[2]:
            symbols = str(records[2]).split(',')
            oryGene_ds[oryGeneID]["Symbols"].extend(symbols)
        if records[4]:
            alt_names = str(records[4]).split(',')
            oryGene_ds[oryGeneID]["Alt_names"].extend(alt_names)
        if records[9]:
            traits = str(records[9]).split(',')
            oryGene_ds[oryGeneID]["Trait_class"].extend(traits)
        if records[10]:
            # re.split(), not str.split(), on a multi-character pattern --
            # str.split() treats ',|\|/' as a literal substring (it is not a
            # regex), so a multi-valued field like "Os06g0157800,Os06g0158100"
            # never actually got split: the whole string stayed one token,
            # never matched rap_pattern below, and both RAP ids were silently
            # dropped. Verified against a real export: 1,276 genes affected,
            # 1,455 RAP id links lost. MSU_id/Gramene_id below already used
            # re.split() correctly; this just brings RAP_id in line with them.
            rapIds = re.split(r',|\||\/', str(records[10]))
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
            oryGene_ds[oryGeneID]["MSU_id"].extend(ids)
        if records[12]:
            records[12] = re.sub('\\r$', '', str(records[12]))
            gr_ids = re.split(r'/|,', str(records[12]))
            gr_ids = [x.rstrip() for x in gr_ids]
            oryGene_ds[oryGeneID]["Gramene_id"].extend(gr_ids)
    return oryGene_ds


def oryzaBaseRDF(infile, output_file, pub_dict, pub_dict_rapdb, pub_dict_msu,
                  pub_dict_hub, dict_po, dict_to):
    print("*************Parsing OryzaBase gene data ***********\n")

    orygene_ds = oryzaBaseParser(infile)
    gene_count = len(orygene_ds)

    print("Number of genes: %s\n" % (str(gene_count)))
    print("OryzaBase gene data has been parsed!\n")
    print("*************************************\n\n")

    ttl_handle = open(output_file, "w")
    pub_handle = open(pub_dict, 'w')
    pub_handle_rapdb = open(pub_dict_rapdb, 'w')
    pub_handle_msu = open(pub_dict_msu, 'w')
    pub_handle_genomeHub = open(pub_dict_hub, 'w')
    PO_handle_genomeHub = open(dict_po, 'w')
    TO_handle_genomeHub = open(dict_to, 'w')

    print("************* OryzaBase RDF conversion begins***********\n")

    ttl_handle.write(str(getRDFHeaders()))
    for oryid in orygene_ds:
        rap_id = ''
        ttl_buffer = ''
        pub_buffer = ''
        ttl_buffer += oryzabase_ns + oryid + "\n"
        ttl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
        ttl_buffer += "\t" + dcterms_ns + "identifier" + "\t" + '"%s"' % turtle_literal(oryid) + " ;\n"
        for item in orygene_ds[oryid]:
            if item == 'Reco_symbol':
                if orygene_ds[oryid][item]:
                    label = re.sub('^\s+|^\t+|\s+$', '', orygene_ds[oryid][item])
                    ttl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % turtle_literal(label) + " ;\n"
                    ttl_buffer += "\t" + skos_ns + "prefLabel" + "\t" + '"%s"' % turtle_literal(label) + " ;\n"
                    if label != '_':
                        pub_buffer += label + "\t" + oryzabase_uri + oryid + "\n"
                    pubRAPDB_writer(label, orygene_ds, oryid, pub_handle_rapdb)
                    pubMSU_writer(label, orygene_ds, oryid, pub_handle_msu)
            if item == 'Reco_name':
                if orygene_ds[oryid][item]:
                    description = orygene_ds[oryid][item]
                    ttl_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % turtle_literal(description) + " ;\n"
            if item == 'Name':
                if orygene_ds[oryid][item]:
                    altlabel = orygene_ds[oryid][item]
                    ttl_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % turtle_literal(altlabel) + " ;\n"
            if item == 'Explanation':
                if orygene_ds[oryid][item]:
                    comment = orygene_ds[oryid][item]
                    ttl_buffer += "\t" + rdfs_ns + "comment" + "\t" + '"%s"' % turtle_literal(comment) + " ;\n"
            if item == 'Chromosome':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "chromosome" + "\t" + '"%s"' % turtle_literal(orygene_ds[oryid][item]) + " ;\n"
            if item == 'Mutant':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasMutant" + "\t" + '"%s"' % turtle_literal(orygene_ds[oryid][item]) + " ;\n"
            if item == 'Arm':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasChromosomeArm" + "\t" + '"%s"' % turtle_literal(orygene_ds[oryid][item]) + " ;\n"
            if item == 'Locus':
                if orygene_ds[oryid][item]:
                    ttl_buffer += "\t" + base_vocab_ns + "hasLocus" + "\t" + '"%s"' % turtle_literal(orygene_ds[oryid][item]) + " ;\n"
            if item == 'Symbols':
                if orygene_ds[oryid][item]:
                    for symbol in orygene_ds[oryid][item]:
                        if symbol is not None:
                            symbol = re.sub(';', '/', symbol)
                            symbol = re.sub('^\s+|^\t+|\n+|\s+$', '', symbol)
                            if symbol != '_':
                                ttl_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % turtle_literal(symbol) + " ;\n"
                                pub_buffer += symbol + "\t" + oryzabase_uri + oryid + "\n"
                                pubRAPDB_writer(symbol, orygene_ds, oryid, pub_handle_rapdb)
                                pubMSU_writer(symbol, orygene_ds, oryid, pub_handle_msu)
            if item == 'Alt_names':
                if orygene_ds[oryid][item]:
                    for alt_name in orygene_ds[oryid][item]:
                        alt_name = re.sub('^\s+|^\t+|\n+|\s+$', '', alt_name)
                        if alt_name != '_':
                            ttl_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % turtle_literal(alt_name) + " ;\n"
            if item == 'RAP_id':
                if orygene_ds[oryid][item]:
                    for rap_id in orygene_ds[oryid][item]:
                        if re.match(rap_pattern, rap_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + re.sub('\s+', '', rap_id) + " ;\n"
                            ttl_buffer += "\t" + owl_ns + "sameAs" + "\t" + res_ns + re.sub('\s+', '', rap_id) + " ;\n"
                            pub_buffer += rap_id + "\t" + rapdb_gene_uri + re.sub('\s+', '', rap_id) + "\n"
                            pub_handle_genomeHub.write(re.sub('\s+', '', rap_id) + "\t" + ''.join(orygene_ds[oryid]['Reco_symbol']) + "\t" + ''.join(orygene_ds[oryid]['Symbols']) + "\n")
            if item == 'MSU_id':
                if orygene_ds[oryid][item]:
                    for msu_id in orygene_ds[oryid][item]:
                        if re.match(tigr_pattern, msu_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + res_ns + re.sub('\s+', '', msu_id) + " ;\n"
                            ttl_buffer += "\t" + owl_ns + "sameAs" + "\t" + res_ns + re.sub('\s+', '', msu_id) + " ;\n"
                            pub_buffer += msu_id + "\t" + sniplay_gene_uri + re.sub('\.\d\s*$', '', msu_id) + "\n"
            if item == 'Gramene_id':
                if orygene_ds[oryid][item]:
                    for gr_id in orygene_ds[oryid][item]:
                        if re.match(gramene_pattern, gr_id):
                            ttl_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + gr_g_ns + re.sub('\s+', '', gr_id) + " ;\n"
            if item == 'Protein_name':
                if orygene_ds[oryid][item]:
                    if re.match(prot_pattern, orygene_ds[oryid][item]):
                        ttl_buffer += "\t" + sio_ns + "SIO_010078" + "\t" + '"%s"' % turtle_literal(orygene_ds[oryid][item]) + " ;\n"  # encodes
            if item == 'Trait_class':
                if orygene_ds[oryid][item]:
                    for traits in orygene_ds[oryid][item]:
                        for trait in traits.split('-'):
                            trait = re.sub('^\s+|\s+$', '', trait)
                            ttl_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + '"%s"' % turtle_literal(trait) + " ;\n"
            if item == 'TO_id':
                if orygene_ds[oryid][item]:
                    for to_term in orygene_ds[oryid][item]:
                        for to_id in to_term.split('-'):
                            to_id = re.sub(':', '_', to_id)
                            to_id = re.sub('^\s|\s$', '', to_id)
                            if re.match(r"^TO_[0-9]{7}$", to_id):
                                ttl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + obo_ns + to_id + " ;\n"
                        if rap_id:
                            TO_handle_genomeHub.write(rap_id + "\t" + re.sub('^\s|\s$', '', to_term) + "\n")
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
                                ttl_buffer += "\t" + obo_ns + "RO_0002206" + "\t" + obo_ns + po_id + " ;\n"  # expressed in
                        if rap_id:
                            PO_handle_genomeHub.write(rap_id + "\t" + re.sub('^\s|\s$', '', po_term) + "\n")

        ttl_buffer = re.sub(' ;$', ' .\n', ttl_buffer)
        pub_handle.write(pub_buffer)
        ttl_handle.write(ttl_buffer)

    ttl_handle.close()
    print("************* OryzaBase RDF completed ************!\n\n")


def pubRAPDB_writer(label, orygene_ds, oryid, pub_handle_rapdb):
    if orygene_ds[oryid]['RAP_id']:
        for rap_id in orygene_ds[oryid]['RAP_id']:
            rap_id = re.sub('^\s+|^\t+|\s+$', '', rap_id)
            if rap_id.split('/'):
                rap_dic = rap_id.split('/')
                for rap_id in rap_dic:
                    if re.match(rap_pattern, rap_id):
                        if label != '_':
                            pub_handle_rapdb.write(label + "\t" + rapdb_gene_uri + rap_id + "\n")
            else:
                if re.match(rap_pattern, rap_id):
                    if label != '_':
                        pub_handle_rapdb.write(label + "\t" + rapdb_gene_uri + rap_id + "\n")


def pubMSU_writer(label, orygene_ds, oryid, pub_handle_msu):
    if orygene_ds[oryid]['MSU_id']:
        for msu_id in orygene_ds[oryid]['MSU_id']:
            msu_id = re.sub('^\s+|^\t+|\s+$', '', msu_id)
            msu_id = re.sub('\.\d\s*$', '', msu_id)
            if re.match(tigr_pattern, msu_id):
                if label != '_':
                    pub_handle_msu.write(label + "\t" + sniplay_gene_uri + msu_id + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert OryzaBase's gene list export to RDF, plus the "
                    "side mapping files (RAP-DB/MSU publication dictionaries "
                    "and the genome hub's PO/TO/gene-symbol tables) other "
                    "AgroLD_ETL parsers join against.",
        epilog="Example:\n"
               "  %(prog)s OryzabaseGeneListEn_20261022010113.txt oryzabase_genes.ttl",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_file",
                        help="OryzaBase gene list export, e.g. "
                             "OryzabaseGeneListEn_YYYYMMDDHHMMSS.txt")
    parser.add_argument("output_file", help="output Turtle file for the gene RDF")
    parser.add_argument("--pub-dict", default=None,
                        help="symbol/RAP-db-id publication dictionary "
                             "(default: pub_dictionnary.txt next to output_file)")
    parser.add_argument("--pub-dict-rapdb", default=None,
                        help="default: pub_dictionnary_rapdb.txt next to output_file")
    parser.add_argument("--pub-dict-msu", default=None,
                        help="default: pub_dictionnary_msu.txt next to output_file")
    parser.add_argument("--pub-dict-hub", default=None,
                        help="default: pub_genome_hub.txt next to output_file")
    parser.add_argument("--dict-po", default=None,
                        help="default: dict_PO_hub.txt next to output_file")
    parser.add_argument("--dict-to", default=None,
                        help="default: dict_TO_hub.txt next to output_file")
    args = parser.parse_args()

    output_dir = os.path.dirname(args.output_file) or "."
    pub_dict = args.pub_dict or os.path.join(output_dir, "pub_dictionnary.txt")
    pub_dict_rapdb = args.pub_dict_rapdb or os.path.join(output_dir, "pub_dictionnary_rapdb.txt")
    pub_dict_msu = args.pub_dict_msu or os.path.join(output_dir, "pub_dictionnary_msu.txt")
    pub_dict_hub = args.pub_dict_hub or os.path.join(output_dir, "pub_genome_hub.txt")
    dict_po = args.dict_po or os.path.join(output_dir, "dict_PO_hub.txt")
    dict_to = args.dict_to or os.path.join(output_dir, "dict_TO_hub.txt")

    oryzaBaseRDF(args.input_file, args.output_file, pub_dict, pub_dict_rapdb,
                 pub_dict_msu, pub_dict_hub, dict_po, dict_to)


if __name__ == "__main__":
    main()
