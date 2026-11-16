import argparse
import datetime
import gzip
import os
import re
import sys
import yaml
from urllib.parse import urljoin
from urllib.request import Request, urlopen

TF_LIST_PAGE = "https://planttfdb.gao-lab.org/download.php"
TF_LIST_BASE_URL = "https://planttfdb.gao-lab.org/"

# The TF-target regulatory network isn't on PlantTFDB itself -- it's on its
# sister site PlantRegMap, section "3. Regulations/Networks".
REGULATION_PAGE = "https://plantregmap.gao-lab.org/download.php"
REGULATION_BASE_URL = "https://plantregmap.gao-lab.org/"

# Matches one row of the TF list table on TF_LIST_PAGE, e.g.:
#   <tr><td><a href="/index.php?sp=Ath"><i>Arabidopsis thaliana</i></a></td>
#   <td><a href="download/TF_list/Ath_TF_list.txt.gz">download</a></td>...
TF_LIST_ROW = re.compile(
    r'<tr><td><a href="/index\.php\?sp=[A-Za-z]+"><i>([^<]+)</i></a></td>'
    r'<td><a href="(download/TF_list/[^"]+)">')

# Matches one row of the "Regulations/Networks" table on REGULATION_PAGE,
# e.g.:
#   <tr><td class="left"><i>Arabidopsis thaliana</i></td>
#   <td class="left"><a href="./download_ftp.php?filepath=...regulation_from_motif_Ath.txt">download</a></td>
#   <td class="left">-</td><td class="left">-</td>
#   <td class="left"><a href="./download_ftp.php?filepath=...regulation_merged_Ath.txt">download</a></td></tr>
# Only the last ("Regulation merged") column's link is wanted; the other
# three columns are sometimes "-" when that breakdown isn't available.
REGULATION_ROW = re.compile(
    r'<tr><td class="left"><i>([^<]+)</i></td>(?P<rest>(?:(?!</tr>).)*?)</tr>',
    re.DOTALL)
REGULATION_MERGED_HREF = re.compile(r'href="(\./download_ftp\.php\?filepath=[^"]*regulation_merged_[^"]+)"')

DEFAULT_SPECIES_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "ensembl_rdf", "ensembl_rdf", "config", "species_agrold.yaml")


def log(msg):
    print(f'[{datetime.datetime.now()}] {msg}', file=sys.stderr)


def fetch(url):
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as res:
        return res.read()


def load_species_list(species_file):
    with open(species_file) as f:
        return yaml.safe_load(f)["species"]


def production_name_to_latin(name):
    """arabidopsis_thaliana -> Arabidopsis thaliana

    Only handles the common two-word case; production names with an assembly
    suffix (e.g. gossypium_raimondii_gca025698545v1rs) won't match anything
    on PlantTFDB and are reported as skipped, not an error -- PlantTFDB
    doesn't cover every AgroLD species.
    """
    genus, _, epithet = name.partition("_")
    return f"{genus.capitalize()} {epithet}"


def parse_tf_list_table(html):
    """{latin name: download/TF_list/... href} from PlantTFDB's download page."""
    text = html.decode("utf-8", errors="replace")
    mapping = {latin: href for latin, href in TF_LIST_ROW.findall(text)}
    if not mapping:
        raise RuntimeError(
            "No species rows found in the TF list table -- PlantTFDB's "
            "download page layout may have changed")
    return mapping


def parse_regulation_table(html):
    """{latin name: ./download_ftp.php?filepath=...regulation_merged_...}
    from PlantRegMap's "Regulations/Networks" table. A species with no
    merged-regulation column (just "-") is left out of the mapping."""
    text = html.decode("utf-8", errors="replace")
    mapping = {}
    for latin, rest in REGULATION_ROW.findall(text):
        m = REGULATION_MERGED_HREF.search(rest)
        if m:
            mapping[latin] = m.group(1)
    if not mapping:
        raise RuntimeError(
            "No species rows found in the Regulations/Networks table -- "
            "PlantRegMap's download page layout may have changed")
    return mapping


def download_file(url, dest, gzipped):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    log(f'Downloading {url}')
    data = fetch(url)
    if gzipped:
        import io
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
            data = gz.read()
    with open(dest + ".part", "wb") as f:
        f.write(data)
    os.replace(dest + ".part", dest)
    log(f'Wrote {dest}')


def main():
    parser = argparse.ArgumentParser(
        description="Download PlantTFDB's TF list and PlantRegMap's merged "
                    "TF-target regulation network for the species AgroLD "
                    "covers, one subdirectory per species under "
                    "--output-dir, named and laid out as config.yaml's "
                    "planttfdb section and the Snakefile expect.",
        epilog="Example:\n"
               "  %(prog)s -s arabidopsis_thaliana",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--species-file", default=DEFAULT_SPECIES_FILE,
                        help="YAML file with a `species` list (default: "
                             "ensembl_rdf's species_agrold.yaml, the same "
                             "file the Snakefile resolves taxon ids from; "
                             "assumed to be a sibling checkout)")
    parser.add_argument("-s", "--species", nargs="+", metavar="NAME",
                        help="only download these species (production name, "
                             "e.g. arabidopsis_thaliana); default: every "
                             "species in --species-file that PlantTFDB has")
    parser.add_argument("-o", "--output-dir", default="data/planttfdb",
                        help="base directory to write per-species "
                             "subdirectories into (default: data/planttfdb)")
    parser.add_argument("--tf-list-filename", default="Mes_TF_list.txt",
                        help="name to give the downloaded TF list in each "
                             "species subdirectory (default: Mes_TF_list.txt, "
                             "matching config.yaml's planttfdb.tf_list)")
    parser.add_argument("--regulation-filename", default="regulation_merged_Mes.txt",
                        help="name to give the downloaded regulation network "
                             "in each species subdirectory (default: "
                             "regulation_merged_Mes.txt, matching "
                             "config.yaml's planttfdb.regulation_file)")
    parser.add_argument("--skip-tf-list", action="store_true",
                        help="don't download the TF list")
    parser.add_argument("--skip-regulation", action="store_true",
                        help="don't download the regulation network "
                             "(it's much larger than the TF list -- tens "
                             "of MB per species)")
    args = parser.parse_args()

    wanted = args.species or load_species_list(args.species_file)
    log(f'{len(wanted)} species requested')

    tf_list_by_latin = {} if args.skip_tf_list else parse_tf_list_table(fetch(TF_LIST_PAGE))
    regulation_by_latin = {} if args.skip_regulation else parse_regulation_table(fetch(REGULATION_PAGE))

    downloaded, skipped = [], []
    for name in wanted:
        latin = production_name_to_latin(name)
        got_something = False

        if not args.skip_tf_list:
            href = tf_list_by_latin.get(latin)
            if href is not None:
                download_file(TF_LIST_BASE_URL + href,
                               os.path.join(args.output_dir, name, args.tf_list_filename),
                               gzipped=True)
                got_something = True

        if not args.skip_regulation:
            href = regulation_by_latin.get(latin)
            if href is not None:
                download_file(urljoin(REGULATION_BASE_URL, href),
                               os.path.join(args.output_dir, name, args.regulation_filename),
                               gzipped=False)
                got_something = True

        (downloaded if got_something else skipped).append(name)

    log(f'Downloaded {len(downloaded)} species: {downloaded}')
    if skipped:
        log(f'Skipped (no PlantTFDB/PlantRegMap entry for): {len(skipped)} species: {skipped}')

    if not downloaded:
        # Nothing succeeded -- in particular the case Snakemake's
        # download_planttfdb rule relies on, where -s names exactly one
        # species: exit non-zero so a species with no PlantTFDB/PlantRegMap
        # coverage shows up as a clear job failure instead of a silently-
        # missing output file.
        sys.exit(1)


if __name__ == "__main__":
    main()
