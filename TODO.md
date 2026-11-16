# TODO

## Infer gene --encodes--> protein via SPARQL CONSTRUCT (UniProt OrderedLocusNames x Ensembl gene identifier)

### Context

`mappings_<taxon>_<bank>_gene_prot.tsv` (a side output of `riceKB/uniprotToRDF.py`) maps a UniProt
accession to its `OrderedLocusNames` symbol(s), but nothing currently converts this into RDF.

Checked empirically against real Ensembl core dumps (rice taxon 39947, Arabidopsis taxon 3702,
`ensembl_rdf_data/{oryza_sativa,arabidopsis_thaliana}_core_*`):
- Ensembl's own xref already links far more UniProt accessions than this mapping does (e.g. rice:
  58,080 via xref vs 4,147 in the TSV), but only as a generic `rdfs:seeAlso`.
- ~95-97% of the TSV's accessions are already covered by that xref, but the mapping would still add
  a semantically precise `agrold_vocabulary:encodes` relation (already used elsewhere, e.g.
  `riceKB/grameneParsersAPI.py`) in place of the generic `seeAlso`, plus a genuinely new ~5% of
  accessions/species where Ensembl's xref has nothing at all (780 for rice, 555 for Arabidopsis).
- The UniProt `OrderedLocusNames` symbol is (case-insensitively) the real Ensembl gene `stable_id`.

### Decision

Instead of a new Python script joining against a TSV, do this with SPARQL directly on the loaded
store: both literals already exist in the data once loaded --
- `ensembl_rdf` already emits `<gene> dcterms:identifier "<stable_id>"` on every gene
  (`ensembl_rdf/rdf_converter_ensembl_db.py:349`)
- `riceKB/uniprotToRDF.py` already emits `<protein> agrold_vocabulary:orderedLocusName "<symbol>"`
  on every protein (fixed below)

A `SPARQL CONSTRUCT ... INSERT` run once both datasets are loaded, joining the two literals
case-insensitively, materializes:
```
?gene agrold_vocabulary:encodes ?protein .
```
(Considered a SHACL-AF rule (`sh:SPARQLRule`) for this instead of a plain CONSTRUCT -- not pursuing
it: SHACL-AF ("rules", as opposed to SHACL Core validation) is an optional extension, and it's not
confirmed the target Virtuoso instance executes/materializes SHACL rules rather than only validating
conformance. A plain SPARQL CONSTRUCT/INSERT needs no such support and does the same job.)

### Done already

- [x] `riceKB/uniprotToRDF.py`: UniProt sometimes merges two loci into one `OrderedLocusNames` entry
  with a `/` (e.g. `Os02g0318450/Os02g0318500`) instead of two comma-separated entries. Fixed to
  split on `/` too, emitting one `altSymbol`/`orderedLocusName` literal (and one `gene_prot.tsv` row)
  per actual locus, so a literal join matches an individual locus id instead of the combined string.

### To do

- [ ] Write the SPARQL CONSTRUCT/INSERT (case-insensitive join on `dcterms:identifier` /
  `agrold_vocabulary:orderedLocusName`, predicate `agrold_vocabulary:encodes`) and decide where it
  runs from (one-off `isql`/SPARQL-endpoint step after a full load, likely alongside
  `bin/prepare_virtuoso_import.py`'s `load.sql`, not a new Snakemake rule since it needs both
  datasets already loaded together).
