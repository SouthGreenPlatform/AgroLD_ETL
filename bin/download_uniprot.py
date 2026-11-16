import argparse
import datetime
import os
import sys
from urllib.request import Request, urlopen

# ftp.uniprot.org itself didn't respond from where this was developed;
# EBI's mirror carries the exact same taxonomic_divisions layout and files.
UNIPROT_FTP_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
    "knowledgebase/taxonomic_divisions/")


CHUNK_SIZE = 1024 * 1024  # 1MB -- keeps memory use flat regardless of file size


def log(msg):
    print(f'[{datetime.datetime.now()}] {msg}', file=sys.stderr)


def download(url, dest):
    """Stream url straight to dest, still gzip-compressed -- riceKB/
    uniprotToRDF.py decompresses on the fly as it reads, so there is no
    need (and, for uniprot_trembl_plants.dat.gz's ~14GB, no room) to ever
    hold the whole response or the decompressed content in memory or on
    disk at once.
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    log(f'Downloading {url}')
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as res, open(dest + ".part", "wb") as f:
        while chunk := res.read(CHUNK_SIZE):
            f.write(chunk)
    os.replace(dest + ".part", dest)
    log(f'Wrote {dest}')


def main():
    parser = argparse.ArgumentParser(
        description="Download a UniProt taxonomic division flat-file (e.g. "
                    "uniprot_sprot_plants.dat.gz) from UniProt's FTP site. "
                    "This is ONE file shared by every species in that "
                    "division -- unlike STRING/PlantTFDB, there is no "
                    "per-species variant to pick, so riceKB/uniprotToRDF.py "
                    "filters the species it wants out of this file itself.",
        epilog="Example:\n"
               "  %(prog)s --bank sprot --division plants",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bank", default="sprot", choices=["sprot", "trembl"],
                        help="sprot (Swiss-Prot, reviewed) or trembl "
                             "(unreviewed, much larger per division) "
                             "(default: sprot)")
    parser.add_argument("--division", default="plants",
                        help="UniProt taxonomic division name, as used in "
                             "the UniProt FTP filenames (default: plants)")
    parser.add_argument("-o", "--output-dir", default="data/uniprot",
                        help="directory to write the still-gzipped .dat.gz "
                             "file into (default: data/uniprot)")
    args = parser.parse_args()

    filename = f"uniprot_{args.bank}_{args.division}.dat.gz"
    url = UNIPROT_FTP_BASE + filename
    dest = os.path.join(args.output_dir, filename)
    download(url, dest)


if __name__ == "__main__":
    main()
