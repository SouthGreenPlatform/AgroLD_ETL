#!/usr/bin/env python
'''
Created on May, 2017
Updated on Dec, 2021
The qtaroModel module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Qtaro data

Input:
    Gene export ('|'-separated, header row, verified against a real export):
        ogro_id|gene|gene_symbol|character_major|character_minor|chromosome|start|end|locus_id|browse|isolation|objective|doi
    QTL export ('|'-separated, header row):
        id|qtl_gene_name|character_major|character_minor|marker_volume|chromosome|start|end|source|lod|crossed_a|crossed_b|direction|references|ref_number|marker|marker_physical|marker_fine_1|marker_fine_2|marker_fine_3|marker_interval_1|marker_interval_2|marker_interval_3|marker_single|analytic_group|group_size|trait_name|contribution_rate|synergy_effect|publication_year

TODO:
    1) Add documentation
    2) better Error handling
@author: larmande
'''
import argparse
import re
import numpy as np
import pandas as pd
from riceKB.globalVars import *
from riceKB.utils import *

taxon_id = "39947"


def geneParser(infile):
    """Qtaro gene export -> DataFrame, dropping rows with no locus_id (RAP
    id) -- those can't be resolved to an AgroLD gene resource.
    """
    array = pd.read_csv(infile, sep="|", dtype='str', quotechar='"', encoding="ISO-8859-1")
    array['locus_id'] = array['locus_id'].replace('', np.nan)
    return array.dropna(subset=['locus_id'])


def qtlParser(infile):
    """Qtaro QTL export -> DataFrame (accessed positionally, see module docstring).

    Unlike the gene export, this file has no header row (verified against a
    real export, datasets/qtaro.csv) -- header=None keeps pandas from
    treating the first data row as column names.
    """
    return pd.read_csv(infile, sep="|", dtype='str', header=None, encoding="ISO-8859-1")


'''
 RDF Converters
'''
def qtaroGeneRDF(infile, output_file, msu_rap_mapping=None):
    gene_counter = 0
    converted_counter = 0
    unmapped_counter = 0
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')

    print("*********** Parsing Qtaro Gene data ***************\n")
    gene_ds = geneParser(infile)

    print("************* Qtaro Gene RDF conversion begins***********\n")
    with open(output_file, "w") as outHandle:
        outHandle.write(str(getRDFHeaders()))
        for records in gene_ds.to_numpy():
            locus_id = records[8]
            if not rap_pattern.match(locus_id):
                # A handful of QTARO rows give an MSU id (e.g.
                # "LOC_Os03g15680") instead of a RAP/IRGSP one -- convert it
                # via RAP-DB's RAP-MSU mapping rather than dropping the gene.
                rap_id = msu_rap_mapping.get(locus_id.split('.')[0]) if msu_rap_mapping else None
                if rap_id is None:
                    unmapped_counter += 1
                    continue
                locus_id = rap_id
                converted_counter += 1
            gene_counter += 1

            # The source wraps a field's whole value in an extra pair of
            # quotes (-> '""...""' in the raw file, unescaped by the CSV
            # parser into a literal leading/trailing '"') whenever the text
            # itself contains a character the exporter treats as sensitive
            # (seen on gene/objective, and once on doi together with a
            # comma-separated second DOI) -- strip that wrapper before
            # escaping the value for Turtle, or the literal quote survives
            # into the output and breaks the syntax.
            gene_buffer = base_resource_ns + locus_id + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % turtle_literal(records[1].strip('"')) + " ;\n"
            gene_buffer += "\t" + skos_ns + "prefSymbol" + "\t" + '"%s"' % turtle_literal(records[2]) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % turtle_literal(records[3]) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % turtle_literal(records[4]) + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + ensembl_ns + locus_id + " ;\n"
            gene_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % turtle_literal(records[11].strip('"')) + " ;\n"
            if pd.notna(records[12]):
                for doi in records[12].strip('"').split(','):
                    doi = doi.strip()
                    if doi:
                        gene_buffer += "\t" + dc_ns + "references" + "\t" + "<" + doi_uri + doi + "> ;\n"
            gene_buffer = re.sub(' ;$', ' .\n', gene_buffer)
            outHandle.write(gene_buffer)
    print("Number of Genes: %s (%s converted from an MSU id, %s dropped for "
          "lacking any RAP id and having no MSU match)\n"
          % (gene_counter, converted_counter, unmapped_counter))
    print("********* Qtaro GENE RDF completed ***********\n")


def qtaroQTLRDF(infile, output_file):
    qtl_counter = 0

    print("*********** Parsing Qtaro QTL data ***************\n")
    qtl_ds = qtlParser(infile)

    print("************* Qtaro QTL RDF conversion begins***********\n")
    with open(output_file, "w") as outHandle:
        outHandle.write(str(getRDFHeaders()))
        for records in qtl_ds.to_numpy():
            qtl_counter += 1

            # Same quote-wrapping convention as the gene export (see
            # qtaroGeneRDF) shows up here on qtl_gene_name and trait_name --
            # strip it before escaping, or the literal quote breaks the
            # Turtle syntax.
            qtl_buffer = "<" + base_resource_uri + "qtaro.qtl/" + records[0] + ">\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "QTL" + " ;\n"
            if pd.notna(records[1]):
                qtl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % turtle_literal(records[1].strip('"')) + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "location" + "\t" + "<" + chromosome_uri + taxon_id + "/{}:{}-{}:1>".format(records[5], records[6], records[7]) + " ;\n"
            if pd.notna(records[2]):
                qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % turtle_literal(records[2]) + " ;\n"
            if pd.notna(records[3]):
                qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % turtle_literal(records[3]) + " ;\n"
            if pd.notna(records[26]):
                qtl_buffer += "\t" + base_vocab_ns + "hasTrait" + "\t" + '"%s"' % turtle_literal(records[26].strip('"')) + " ;\n"
            if pd.notna(records[9]):
                qtl_buffer += "\t" + base_vocab_ns + "lod" + "\t" + '"%s"' % turtle_literal(records[9]) + " ;\n"
            qtl_buffer += "\t" + base_vocab_ns + "date" + "\t" + '"%s"' % turtle_literal(records[29]) + " .\n\n"

            # Region
            qtl_buffer += "<" + chromosome_uri + taxon_id + "/{}:{}-{}:1>".format(records[5], records[6], records[7]) + "\n"
            qtl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % turtle_literal(taxon_id + "/{}:{}-{}:1".format(records[5], records[6], records[7])) + " ;\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "begin" + "\t" + "<" + chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5], records[6]) + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "end" + "\t" + "<" + chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5], records[7]) + " .\n\n"

            # Position 1
            qtl_buffer += "<" + chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5], records[6]) + "\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ForwardStrandPosition" + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "position" + "\t" + records[6] + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + "/{}>".format(records[5]) + " .\n\n"

            # Position 2
            qtl_buffer += "<" + chromosome_uri + taxon_id + "/{}:{}:1>".format(records[5], records[7]) + "\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
            qtl_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ForwardStrandPosition" + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "position" + "\t" + records[7] + " ;\n"
            qtl_buffer += "\t" + faldo_ns + "reference" + "\t" + "<" + chromosome_uri + taxon_id + "/{}>".format(records[5]) + " .\n\n"

            outHandle.write(qtl_buffer)
    print("Number of QTLs: %s\n" % (str(qtl_counter)))
    print("********* Qtaro QTL RDF completed ***********\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a QTARO gene or QTL export to RDF.")
    sub = parser.add_subparsers(dest="command", required=True)

    gene = sub.add_parser("genes", help="convert the QTARO gene export")
    gene.add_argument("input_file", help="QTARO gene export ('|'-separated, header row)")
    gene.add_argument("output_file", help="output Turtle file")
    gene.add_argument("--msu-rap-mapping",
                       help="RAP-DB's RAP-MSU mapping file (e.g. RAP-MSU_2026-02-05.txt.gz, "
                            "gzip-compressed or plain text), "
                            "used to convert the handful of rows that give an MSU id "
                            "(LOC_Os...) instead of a RAP/IRGSP one; without it, those "
                            "rows are dropped instead of converted")

    qtl = sub.add_parser("qtl", help="convert the QTARO QTL export")
    qtl.add_argument("input_file", help="QTARO QTL export ('|'-separated, header row)")
    qtl.add_argument("output_file", help="output Turtle file")

    args = parser.parse_args()
    if args.command == "genes":
        mapping = load_msu_rap_mapping(args.msu_rap_mapping) if args.msu_rap_mapping else None
        qtaroGeneRDF(args.input_file, args.output_file, mapping)
    else:
        qtaroQTLRDF(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
