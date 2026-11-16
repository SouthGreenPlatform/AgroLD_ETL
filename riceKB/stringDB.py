#!/usr/bin/env python
'''
Created on June 8, 2020
The stringDB module is created as part of the AgroLD project.

This module contains Parsers, RDF converters and generic functions for handling StringDB data

Input (from https://string-db.org/cgi/download, per species taxon id):
    <taxon_id>.protein.aliases.v<version>.txt  (tab-separated, header row)
        #string_protein_id  alias        source
        3702.A0A0A7EPL0     AT1G08910    KEGG_KEGGID_SHORT
    <taxon_id>.protein.links.v<version>.txt    (space-separated, header row)
        protein1          protein2          combined_score
        3702.A0A0A7EPL0   3702.Q8GYL6       231

STRING's own protein ids (e.g. "A0A0A7EPL0", a UniProt accession for
Arabidopsis) don't match the gene stable ids ensembl_rdf's agrold profile
uses as AgroLD resource URIs, so they can't be used directly to join with
the rest of AgroLD. The aliases file's "KEGG_KEGGID_SHORT" rows give the
native locus id (e.g. "AT1G08910") instead, one per protein id, which for
Ensembl-Plants-derived species is the same string as the Ensembl gene
stable id -- that's what RDFConverter resolves each interacting protein to.
A STRING protein with no KEGG_KEGGID_SHORT alias (e.g. organelle-encoded
genes, which Ensembl Plants doesn't model) can't be resolved this way, and
any links row naming it is skipped (and counted) rather than emitted with
an unresolvable gene id.

#TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: venkatesan

'''
import argparse
import sys
from riceKB.globalVars import *
from riceKB.utils import *


def aliasParser(infile, taxon_id, source="KEGG_KEGGID_SHORT"):
    """{string protein id (taxon prefix stripped): native gene locus id},
    built from the aliases file's rows of the given source. There is at
    most one such row per protein id (verified against real STRING data),
    so the dict stays 1:1 and the last row wins if that ever stops holding.
    """
    prefix = taxon_id + "."
    gene_by_protein = {}
    with open(infile, 'r', encoding='utf-8') as f:
        next(f)  # header: #string_protein_id  alias  source
        for line in f:
            protein_id, alias, alias_source = line.rstrip("\n").split("\t")
            if alias_source != source:
                continue
            if not protein_id.startswith(prefix):
                continue
            gene_by_protein[protein_id[len(prefix):]] = alias
    return gene_by_protein


def RDFConverter(links_file, gene_by_protein, taxon_id, output_file, min_score):
    prefix = taxon_id + "."
    rdf_writer = open(output_file, "w")
    rdf_writer.write(str(getRDFHeaders()))

    line_number = 0
    resolved = 0
    skipped_low_score = 0
    skipped_unresolved = 0
    skipped_wrong_taxon = 0

    print("************* RDF conversion begins***********\n")
    with open(links_file, 'r', encoding='utf-8') as f:
        next(f)  # header: protein1 protein2 combined_score
        for line in f:
            line_number += 1
            protein1, protein2, score = line.split()

            # Checked first and cheaply (int compare, no dict lookup yet) since
            # it drops the bulk of rows -- STRING's own authors call >=800
            # ("high confidence") the threshold worth keeping; below that the
            # network is dominated by low-confidence/textmining-only edges.
            if int(score) < min_score:
                skipped_low_score += 1
                continue

            if not (protein1.startswith(prefix) and protein2.startswith(prefix)):
                skipped_wrong_taxon += 1
                continue

            gene1 = gene_by_protein.get(protein1[len(prefix):])
            gene2 = gene_by_protein.get(protein2[len(prefix):])
            if gene1 is None or gene2 is None:
                skipped_unresolved += 1
                continue

            association = base_resource_uri + "association/string/" + gene1 + "_" + gene2 + "_" + score
            buffer = ''
            buffer += "<" + base_resource_uri + gene1 + ">"
            buffer += "\t" + obo_ns + "RO_0002461" + "\t\t" + "<" + association + "> .\n\n"
            buffer += "<" + association + ">\n"
            buffer += "\t" + rdf_ns + "value" + "\t" + "\"" + score + "\"^^xsd:integer ;\n"
            buffer += "\t" + obo_ns + "RO_0002462" + "\t" + "<" + string_uri + protein1 + "> ;\n"
            buffer += "\t" + obo_ns + "RO_0002463" + "\t" + "<" + string_uri + protein2 + "> ;\n"
            buffer += "\t" + obo_ns + "RO_0002436" + "\t" + "<" + base_resource_uri + gene2 + "> .\n\n"
            rdf_writer.write(buffer)
            resolved += 1

    print("*************** StringDB RDF conversion completed ************\n")
    print(f"{line_number} rows read, {resolved} written, "
          f"{skipped_low_score} skipped (combined_score below --min-score {min_score}), "
          f"{skipped_unresolved} skipped (no gene id for one of the two proteins), "
          f"{skipped_wrong_taxon} skipped (taxon id in file doesn't match --taxon-id)",
          file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Convert a STRING protein-protein association network "
                    "to RDF, mapping STRING's own protein ids to native "
                    "gene locus ids via the aliases file's KEGG_KEGGID_SHORT "
                    "column so the output joins with the rest of AgroLD.")
    parser.add_argument("aliases_file", help="<taxon_id>.protein.aliases.v*.txt")
    parser.add_argument("links_file", help="<taxon_id>.protein.links.v*.txt")
    parser.add_argument("output_file", help="output Turtle file")
    parser.add_argument("--taxon-id", required=True,
                        help="NCBI taxon id of the species this data is for "
                             "(e.g. 3702 for Arabidopsis thaliana) -- also "
                             "the prefix STRING puts on every protein id in "
                             "both input files")
    parser.add_argument("--min-score", type=int, default=800,
                        help="minimum combined_score to keep (0-999, default "
                             "800 -- STRING's own 'high confidence' cutoff); "
                             "lower it to include medium-confidence edges, "
                             "at the cost of a much larger output")
    args = parser.parse_args()

    gene_by_protein = aliasParser(args.aliases_file, args.taxon_id)
    RDFConverter(args.links_file, gene_by_protein, args.taxon_id, args.output_file, args.min_score)


if __name__ == "__main__":
    main()
