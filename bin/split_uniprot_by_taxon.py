#!/usr/bin/env python
"""Split a UniProt taxonomic division flat-file dump (e.g.
uniprot_trembl_plants.dat.gz) into one plain-text file per taxon id, in a
single pass over the input.

riceKB/uniprotToRDF.py filters its one taxon out of the shared division
file by scanning every record in it; run once per taxon (as the Snakefile
used to), that means re-decompressing and re-parsing the whole division
file (up to ~14GB compressed for trembl) once per taxon, which is what
pushed some trembl taxa past even a generous SLURM runtime. This script
does that scan exactly once, fanning each record out to every taxon in
taxon_id that it matches; the uniprot rule then points at its own taxon's
(much smaller) split file instead of the shared dump, with no change to
uniprotToRDF.py itself needed (it already reads plain-text or gzip input
transparently).

Always writes a (possibly empty) output file for every requested taxon id,
even one the input has no record for: riceKB/uniprotToRDF.py's own
exit(1)-on-zero-matches still fires per taxon downstream exactly as it did
before, now reacting to an empty split file instead of to zero matches in
the full dump.
"""
import argparse
import gzip
import os
import re
import sys

NCBI_TAXID_FIELD_RE = re.compile(r'NCBI_TaxID=([\d,\s]+)')


def split(infile, output_dir, taxon_ids):
    """Record boundaries ("//" alone on a line) and each record's taxon
    id(s) (the "OX" line's NCBI_TaxID=...) are read as plain text here,
    matching Bio.SwissProt's own record/taxonomy_id semantics, rather than
    parsing each record with Bio.SwissProt.parse() -- this only needs to
    copy a record's raw text verbatim to the right output file(s), not
    parse its fields just to re-emit the same bytes.
    """
    os.makedirs(output_dir, exist_ok=True)
    wanted = set(taxon_ids)
    handles = {taxon: open(os.path.join(output_dir, f"{taxon}.dat"), "w")
               for taxon in wanted}
    counts = {taxon: 0 for taxon in wanted}

    opener = gzip.open if infile.endswith('.gz') else open
    record_lines = []
    record_taxa = set()
    with opener(infile, 'rt') as f:
        for line in f:
            record_lines.append(line)
            if line.startswith('OX   '):
                match = NCBI_TAXID_FIELD_RE.search(line)
                if match:
                    record_taxa.update(re.findall(r'\d+', match.group(1)))
            if line.rstrip('\n') == '//':
                for taxon in record_taxa & wanted:
                    handles[taxon].writelines(record_lines)
                    counts[taxon] += 1
                record_lines = []
                record_taxa = set()

    for handle in handles.values():
        handle.close()
    for taxon in sorted(wanted, key=int):
        print(f"{taxon}: {counts[taxon]} records", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_file", help="UniProt division dump (.dat or .dat.gz)")
    parser.add_argument("output_dir", help="directory to write <taxon_id>.dat into")
    parser.add_argument("taxon_id", nargs='+', help="NCBI taxon id(s) to split out")
    args = parser.parse_args()
    split(args.input_file, args.output_dir, args.taxon_id)


if __name__ == "__main__":
    main()
