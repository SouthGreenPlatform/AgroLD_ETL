import argparse
import datetime
import gzip
import io
import os
import sys
import urllib.error
import yaml
from urllib.request import Request, urlopen

STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download/"

DEFAULT_SPECIES_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "ensembl_rdf", "ensembl_rdf", "config", "species_agrold.yaml")
DEFAULT_SPECIES_TABLE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "ensembl_rdf", "ensembl_rdf", "config", "species_EnsemblPlants_63.txt")


def log(msg):
    print(f'[{datetime.datetime.now()}] {msg}', file=sys.stderr)


def fetch(url):
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as res:
        return res.read()


def load_agrold_taxon_ids(species_file, species_table):
    """{production_name: taxon_id}, the same lookup the Snakefile's
    load_agrold_taxon_ids() does -- STRING indexes its downloads by NCBI
    taxon id directly, not by species name, so (unlike PlantTFDB/PlantRegMap)
    no name matching or page scraping is needed here at all.
    """
    with open(species_file) as f:
        wanted = set(yaml.safe_load(f)["species"])

    taxon_id = {}
    with open(species_table) as f:
        next(f)  # header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name = fields[1]
            if name in wanted:
                taxon_id[name] = fields[3]

    missing = wanted - set(taxon_id)
    if missing:
        raise RuntimeError(f"{species_table} has no taxon id for: {sorted(missing)}")
    return taxon_id


def download_gzipped(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    log(f'Downloading {url}')
    data = fetch(url)
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
        data = gz.read()
    with open(dest + ".part", "wb") as f:
        f.write(data)
    os.replace(dest + ".part", dest)
    log(f'Wrote {dest}')


def main():
    parser = argparse.ArgumentParser(
        description="Download STRING's protein.aliases and protein.links "
                    "files for the species AgroLD covers, one subdirectory "
                    "per species under --output-dir, named and laid out as "
                    "config.yaml's stringdb section and the Snakefile "
                    "expect.",
        epilog="Example:\n"
               "  %(prog)s -s arabidopsis_thaliana",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--species-file", default=DEFAULT_SPECIES_FILE,
                        help="YAML file with a `species` list (default: "
                             "ensembl_rdf's species_agrold.yaml; assumed to "
                             "be a sibling checkout)")
    parser.add_argument("-t", "--species-table", default=DEFAULT_SPECIES_TABLE,
                        help="tab-separated file mapping each species to its "
                             "NCBI taxon id (default: ensembl_rdf's "
                             "species_EnsemblPlants_63.txt, the same file "
                             "the Snakefile resolves taxon ids from)")
    parser.add_argument("-s", "--species", nargs="+", metavar="NAME",
                        help="only download these species (production name, "
                             "e.g. arabidopsis_thaliana); default: every "
                             "species in --species-file")
    parser.add_argument("-o", "--output-dir", default="data/stringdb",
                        help="base directory to write per-species "
                             "subdirectories into (default: data/stringdb)")
    parser.add_argument("--version", default="12.5",
                        help="STRING release to download (default: 12.5, "
                             "matching config.yaml's stringdb.aliases_suffix "
                             "/ links_suffix -- update both together if you "
                             "change this)")
    args = parser.parse_args()

    taxon_id_by_species = load_agrold_taxon_ids(args.species_file, args.species_table)
    wanted = args.species or sorted(taxon_id_by_species)
    log(f'{len(wanted)} species requested')

    downloaded, skipped = [], []
    for name in wanted:
        taxon_id = taxon_id_by_species.get(name)
        if taxon_id is None:
            log(f'{name}: not in {args.species_file} / {args.species_table}, skipping')
            skipped.append(name)
            continue
        filename_stem = f"{taxon_id}.protein.{{kind}}.v{args.version}.txt"
        try:
            for kind in ("aliases", "links"):
                dataset = f"protein.{kind}.v{args.version}"
                filename = filename_stem.format(kind=kind)
                url = STRING_DOWNLOAD_BASE + dataset + "/" + filename + ".gz"
                dest = os.path.join(args.output_dir, name, filename)
                download_gzipped(url, dest)
            downloaded.append(name)
        except urllib.error.HTTPError as e:
            log(f'{name} (taxon {taxon_id}): {e} -- STRING has no v{args.version} '
                f'data for this taxon id, skipping')
            skipped.append(name)

    log(f'Downloaded {len(downloaded)} species: {downloaded}')
    if skipped:
        log(f'Skipped (no STRING v{args.version} entry for): {len(skipped)} species: {skipped}')

    if not downloaded:
        # Nothing succeeded -- in particular the case Snakemake's
        # download_stringdb rule relies on, where -s names exactly one
        # species: exit non-zero so a species with no STRING coverage shows
        # up as a clear job failure instead of a silently-missing output file.
        sys.exit(1)


if __name__ == "__main__":
    main()
