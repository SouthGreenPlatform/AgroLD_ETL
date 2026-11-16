#!/usr/bin/env python
"""Group every pipeline's Turtle output into one directory per Virtuoso
named graph, ready for Virtuoso's own bulk loader (`ld_dir`) to load it --
one graph per AgroLD species (http://purl.agrold.org/graph/species/<taxon_id>),
one graph per generic (non-species) dataset (http://purl.agrold.org/graph/<name>,
e.g. interpro), and one dedicated graph for Compara orthology
(http://purl.agrold.org/graph/compara/orthology).

Graph directories are populated with symlinks, not copies: several of the
inputs (uniprot's, ensembl_rdf's) are tens of GB, and a copy would double
disk usage for nothing -- `ld_dir` reads every file in a directory (gzipped
or not) into one graph on its own, so there is nothing to concatenate here.

Inputs:
    - AgroLD_ETL's own output_dir (config.yaml): interpro.ttl (generic),
      oryzabase_genes.ttl/qtaro_genes.ttl/qtaro_qtl.ttl (rice-specific,
      taxon RICE_TAXON_ID below), uniprot_<taxon>_<bank>.*.ttl (species, by
      taxon), planttfdb_<species>*.ttl/stringdb_<species>.ttl (species, by
      ensembl_rdf production name), plus any reference ontology (go.owl,
      ncbitaxon.owl, ...) placed in output_dir -- these aren't produced by
      any Snakemake rule, so symlink or copy the downloaded file there by
      hand, like interpro's/oryzabase's manually-placed raw inputs, and add
      it to GENERIC_DATASETS below.
    - ensembl_rdf's sibling checkout: one <species>_core_<version>/
      directory per species (*.ttl/*.ttl.gz), and one *_compara_*/
      orthology.ttl(.gz) -- derived from config.yaml's
      ensembl_rdf.species_file path, or overridable with
      --ensembl-rdf-root.

A *.ttl/*.owl this script doesn't recognize is reported and skipped rather
than guessed at, so a newly-added data source doesn't end up in the wrong
graph (or silently nowhere) until this script is taught about it.

Usage:
    python3 bin/prepare_virtuoso_import.py OUT_DIR
    python3 bin/prepare_virtuoso_import.py --config config.yaml OUT_DIR

Then, on the Virtuoso host (OUT_DIR's directories must be under that
server's DirsAllowed ini setting):
    isql < OUT_DIR/load.sql
"""
import argparse
import os
import re
import sys

# Unlike riceKB/*.py scripts (run by the Snakefile with the repo root on
# PYTHONPATH via shell.prefix), this script is meant to be run by hand, so
# it adds the repo root (this file's parent's parent) to sys.path itself
# instead of relying on the caller's environment.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from riceKB.utils import load_agrold_taxon_ids

GRAPH_BASE = "http://purl.agrold.org/graph"

# OryzaBase and QTARO export gene/QTL data for rice only, but name their
# output generically (oryzabase_genes.ttl, qtaro_genes.ttl, qtaro_qtl.ttl) --
# same taxon riceKB/qtaroModel.py itself already hardcodes, and what the
# shared species table resolves "oryza_sativa" to.
RICE_TAXON_ID = "39947"

RDF_EXTENSIONS = (".ttl", ".owl")

# Reference ontologies (manually placed in output_dir, see this module's
# docstring) go here too, by their exact downloaded filename -- add one
# entry per ontology file rather than matching by extension alone, so an
# unrelated .owl/.ttl dropped into output_dir by mistake is still reported
# as unmatched instead of silently landing in some ontology's graph.
GENERIC_DATASETS = {
    "interpro.ttl": "interpro",
    "go.owl": "go",
    "ncbitaxon.owl": "ncbitaxon",
}

RICE_DATASETS = {"oryzabase_genes.ttl", "reference_oryzabase.ttl", "qtaro_genes.ttl", "qtaro_qtl.ttl"}

UNIPROT_RE = re.compile(r"^uniprot_(?P<taxon>\d+)_(?:sprot|trembl)\..+\.ttl$")
PLANTTFDB_RE = re.compile(r"^planttfdb_(?P<species>.+?)(?:_regulation)?\.ttl$")
STRINGDB_RE = re.compile(r"^stringdb_(?P<species>.+)\.ttl$")
ENSEMBL_RDF_SPECIES_DIR_RE = re.compile(r"^(?P<species>.+)_core_\d.*$")


def species_graph(taxon_id):
    return f"{GRAPH_BASE}/species/{taxon_id}"


def generic_graph(name):
    return f"{GRAPH_BASE}/{name}"


def collect_agrold_etl_files(output_dir, taxon_id_by_species):
    """{graph_iri: [file paths]} from this repo's own output_dir, plus a
    list of *.ttl basenames this function didn't recognize.
    """
    graphs = {}
    unmatched = []

    def add(graph_iri, path):
        graphs.setdefault(graph_iri, []).append(path)

    for name in sorted(os.listdir(output_dir)):
        if not name.endswith(RDF_EXTENSIONS):
            continue  # e.g. mappings_<taxon>_<bank>_gene_prot.tsv -- not RDF
        path = os.path.join(output_dir, name)

        if name in GENERIC_DATASETS:
            add(generic_graph(GENERIC_DATASETS[name]), path)
            continue
        if name in RICE_DATASETS:
            add(species_graph(RICE_TAXON_ID), path)
            continue

        m = UNIPROT_RE.match(name)
        if m:
            add(species_graph(m.group("taxon")), path)
            continue

        m = PLANTTFDB_RE.match(name) or STRINGDB_RE.match(name)
        if m:
            species = m.group("species")
            taxon_id = taxon_id_by_species.get(species)
            if taxon_id is None:
                unmatched.append(name)
                continue
            add(species_graph(taxon_id), path)
            continue

        unmatched.append(name)

    return graphs, unmatched


def collect_compara_files(compara_dir):
    """[file paths] of every orthology.ttl(.gz) directly inside compara_dir
    (e.g. .../ensembl_compara_plants_63_116/orthology.ttl.gz) -- all go into
    the one dedicated compara graph, regardless of where this directory
    lives (a sibling of the species *_core_*/ directories in a plain
    ensembl_rdf checkout, or its own separate path on shared storage, e.g.
    this cluster's /projects/xxl/agrold/compara/<version>/).
    """
    return [os.path.join(compara_dir, name) for name in sorted(os.listdir(compara_dir))
            if name.startswith("orthology.ttl")]


def collect_ensembl_rdf_files(root_dir, taxon_id_by_species):
    """{graph_iri: [file paths]} from the ensembl_rdf sibling checkout:
    every <species>_core_*/'s *.ttl(.gz) into that species' graph, and
    (if present here -- see --compara-dir otherwise) every *_compara_*/'s
    orthology.ttl(.gz) into the one compara graph.
    """
    graphs = {}
    unmatched = []

    def add(graph_iri, path):
        graphs.setdefault(graph_iri, []).append(path)

    for entry in sorted(os.listdir(root_dir)):
        entry_path = os.path.join(root_dir, entry)
        if not os.path.isdir(entry_path):
            continue

        if "_compara_" in entry:
            for path in collect_compara_files(entry_path):
                add(generic_graph("compara/orthology"), path)
            continue

        m = ENSEMBL_RDF_SPECIES_DIR_RE.match(entry)
        if not m:
            continue  # not a *_core_*/*_compara_* directory -- not ours
        species = m.group("species")
        taxon_id = taxon_id_by_species.get(species)
        if taxon_id is None:
            unmatched.append(entry)
            continue
        for name in sorted(os.listdir(entry_path)):
            if name.endswith(".ttl") or name.endswith(".ttl.gz"):
                add(species_graph(taxon_id), os.path.join(entry_path, name))

    return graphs, unmatched


def graph_dir_name(graph_iri):
    """A filesystem-safe directory name for one graph IRI, e.g.
    ".../graph/species/39947" -> "species_39947".
    """
    return graph_iri[len(GRAPH_BASE) + 1:].replace("/", "_")


def populate(graphs, out_dir):
    """One subdirectory of out_dir per graph, symlinking in (not copying)
    every file that belongs to it.
    """
    for graph_iri, paths in sorted(graphs.items()):
        graph_dir = os.path.join(out_dir, graph_dir_name(graph_iri))
        os.makedirs(graph_dir, exist_ok=True)
        for path in paths:
            link = os.path.join(graph_dir, os.path.basename(path))
            target = os.path.abspath(path)
            if os.path.islink(link):
                if os.readlink(link) == target:
                    continue
                raise RuntimeError(
                    f"{link} already links elsewhere ({os.readlink(link)!r} "
                    f"!= {target!r}) -- two different source files share a "
                    "name within the same graph, this script can't tell "
                    "them apart.")
            os.symlink(target, link)


def write_load_script(graphs, out_dir):
    """OUT_DIR/load.sql: one Virtuoso `ld_dir` registration per graph, then
    a single rdf_loader_run() -- see this module's docstring for the
    DirsAllowed caveat before running it on the Virtuoso host.

    The glob is '*', not '*.ttl*': each graph directory is already
    populated with exactly (and only) that graph's files by populate()
    above, including non-.ttl ones like an ontology's .owl, so there's
    nothing left to filter by extension here -- ld_dir sniffs the RDF
    format itself.
    """
    load_sql = os.path.join(out_dir, "load.sql")
    with open(load_sql, "w") as f:
        f.write(
            "-- Generated by bin/prepare_virtuoso_import.py. Each ld_dir's "
            "first argument must be under this Virtuoso server's "
            "DirsAllowed ini setting, or ld_dir itself raises an error.\n\n")
        for graph_iri in sorted(graphs):
            graph_dir = os.path.abspath(os.path.join(out_dir, graph_dir_name(graph_iri)))
            f.write(f"ld_dir('{graph_dir}', '*', '{graph_iri}');\n")
        f.write("\nrdf_loader_run();\ncheckpoint;\n")
    return load_sql


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", help="directory to create the per-graph subdirectories in")
    parser.add_argument("--config", default="config.yaml", help="AgroLD_ETL config (default: config.yaml)")
    parser.add_argument("--ensembl-rdf-root",
                         help="ensembl_rdf sibling checkout root (default: derived from "
                              "config.yaml's ensembl_rdf.species_file); its *_core_*/ "
                              "subdirectories give the species data, and its *_compara_*/ "
                              "one (if present) the compara one, unless --compara-dir is given")
    parser.add_argument("--compara-dir",
                         help="directory holding orthology.ttl(.gz) directly (e.g. "
                              ".../ensembl_compara_plants_63_116), when it isn't a "
                              "*_compara_*/ subdirectory of --ensembl-rdf-root -- e.g. this "
                              "cluster's /projects/xxl/agrold/compara/ensembl_compara_plants_63_116, "
                              "a sibling of /projects/xxl/agrold/ensembl rather than inside it")
    args = parser.parse_args()

    import yaml
    with open(args.config) as f:
        config = yaml.safe_load(f)

    species_file = config["ensembl_rdf"]["species_file"]
    species_table = config["ensembl_rdf"]["species_table"]
    taxon_id_by_species = load_agrold_taxon_ids(species_file, species_table)

    ensembl_rdf_root = args.ensembl_rdf_root or os.path.dirname(
        os.path.dirname(os.path.dirname(species_file)))

    graphs, unmatched_etl = collect_agrold_etl_files(config["output_dir"], taxon_id_by_species)
    rdf_graphs, unmatched_rdf = collect_ensembl_rdf_files(ensembl_rdf_root, taxon_id_by_species)
    for graph_iri, paths in rdf_graphs.items():
        graphs.setdefault(graph_iri, []).extend(paths)

    if args.compara_dir:
        graphs.setdefault(generic_graph("compara/orthology"), []).extend(
            collect_compara_files(args.compara_dir))

    os.makedirs(args.out_dir, exist_ok=True)
    populate(graphs, args.out_dir)
    load_sql = write_load_script(graphs, args.out_dir)

    total_files = sum(len(paths) for paths in graphs.values())
    print(f"{len(graphs)} graphs, {total_files} files linked under {args.out_dir}")
    for graph_iri in sorted(graphs):
        print(f"  {graph_iri}: {len(graphs[graph_iri])} file(s)")
    print(f"Loader script: {load_sql}")

    for name in unmatched_etl:
        print(f"Warning: unrecognized *.ttl in {config['output_dir']}, skipped: {name}", file=sys.stderr)
    for name in unmatched_rdf:
        print(f"Warning: unrecognized entry in {ensembl_rdf_root}, skipped: {name}", file=sys.stderr)


if __name__ == "__main__":
    main()
