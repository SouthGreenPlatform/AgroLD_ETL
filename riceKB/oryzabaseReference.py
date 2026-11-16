#!/usr/bin/env python
"""Link OryzaBase genes to the literature, via OryzaBase's own "Reference"
export -- a different, separately-downloaded export from the gene list
riceKB/oryzaBaseParser.py reads (OryzabaseGeneListEn*.txt).

Only the gene<->PubMed link is emitted (gene dcterms:references pubmed:ID),
the same predicate riceKB/uniprotToRDF.py already uses for the same purpose
on protein resources -- one convention for "this resource is discussed in
PubMed:X" across the pipeline, instead of a second, richer bibo:Article model
that would need its own free-text title/author/journal cleanup and encoding
handling for no reader this pipeline has yet. Rerun with the full
bibliographic fields later if that need shows up.

Input: OryzaBase's "Reference" export (tab-separated, header row), e.g.
Reference_YYYYMMDDHHMMSS.txt, columns (0-indexed, verified against a real
2026-10-22 export -- don't trust an older comment here over this one if the
two ever disagree):
    0  Reference_id
    1  PubMedId
    2  Author
    3  Title
    4  Journal
    5  Volume
    6  Pages
    7  Year
    8  CGSNL Gene Symbol
    9  Gene Name Synonym
Columns 8/9 are comma-separated lists of the gene symbols/synonyms that
publication is about -- the join key against pub_dictionnary.txt (symbol ->
gene URI), itself built by riceKB/oryzaBaseParser.py from the very same
symbols/synonyms in the gene list export. A handful of rows (1 in 28,966 in
that same real export) have stray literal tabs inside a garbled Title field
(pasted MathType markup), which shifts every column after it; such rows are
skipped with a warning rather than guessed at.

@author: larmande
"""
import argparse
import sys
from riceKB.globalVars import *
from riceKB.utils import getRDFHeaders

EXPECTED_FIELDS = 10


def load_pub_dict(pub_dict_file):
    """{symbol: {gene_uri, ...}} from riceKB/oryzaBaseParser.py's
    pub_dictionnary.txt. A set of URIs per symbol because a handful of
    symbols are shared by more than one OryzaBase gene (homonyms) -- link to
    all of them rather than guessing which one a reference meant.
    """
    pub_dict = {}
    with open(pub_dict_file) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                continue
            symbol, gene_uri = fields
            pub_dict.setdefault(symbol, set()).add(gene_uri)
    return pub_dict


def referenceRDF(reference_file, pub_dict_file, output_file):
    pub_dict = load_pub_dict(pub_dict_file)
    links = set()  # (gene_uri, pmid), deduped across rows/symbols/synonyms
    skipped = 0

    with open(reference_file, encoding="utf-8", errors="replace") as f:
        next(f)  # header
        for line_number, line in enumerate(f, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != EXPECTED_FIELDS:
                print(f"Warning: {reference_file}:{line_number} has "
                      f"{len(fields)} fields (expected {EXPECTED_FIELDS}), "
                      "skipped", file=sys.stderr)
                skipped += 1
                continue

            pmid = fields[1].strip()
            if not pmid or not pubmed_pattern.match(pmid):
                continue

            symbols = [s.strip() for s in fields[8].split(",") if s.strip()]
            symbols += [s.strip() for s in fields[9].split(",") if s.strip()]
            for symbol in symbols:
                for gene_uri in pub_dict.get(symbol, ()):
                    links.add((gene_uri, pmid))

    with open(output_file, "w") as out:
        out.write(str(getRDFHeaders()))
        for gene_uri, pmid in sorted(links):
            # gene_uri is a bare full IRI from pub_dictionnary.txt (no
            # angle brackets, no prefix -- it's a plain TSV value, not
            # Turtle), unlike the short oryzabase: prefix riceKB/oryzaBaseParser.py
            # itself uses when writing the same gene as a Turtle subject.
            out.write(f"<{gene_uri}>\n\t{dcterms_ns}references\t{pubmed_ns}{pmid} .\n\n")

    print(f"{len(links)} gene<->PubMed links written to {output_file} "
          f"({skipped} malformed input row(s) skipped)")


def main():
    parser = argparse.ArgumentParser(
        description="Link OryzaBase genes to PubMed via OryzaBase's "
                    "Reference export, joined against the pub_dictionnary.txt "
                    "riceKB/oryzaBaseParser.py already produces.",
        epilog="Example:\n"
               "  %(prog)s Reference_20261022224907.txt pub_dictionnary.txt reference_oryzabase.ttl",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference_file",
                        help="OryzaBase Reference export, e.g. Reference_YYYYMMDDHHMMSS.txt")
    parser.add_argument("pub_dict_file",
                        help="pub_dictionnary.txt, from riceKB/oryzaBaseParser.py's --pub-dict")
    parser.add_argument("output_file", help="output Turtle file")
    args = parser.parse_args()

    referenceRDF(args.reference_file, args.pub_dict_file, args.output_file)


if __name__ == "__main__":
    main()
