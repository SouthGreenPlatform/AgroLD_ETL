"""AgroLD_ETL pipeline, in progress: one rule per data source.

Each source's parser lives in riceKB/ and is invoked as a script (not
imported), so this Snakefile only needs to know its command line, not its
internals. Paths come from config.yaml.

    snakemake --cores 1 --configfile config.yaml
    snakemake --cores 1 -n                       # dry run
    snakemake --profile profiles/slurm           # on the cluster, one sbatch job per rule
    snakemake --profile profiles/slurm uniprot_all  # just uniprot, skipping all's other sources

interpro, oryzaBaseParser, plantTFDB, qtaroModel, stringDB and uniprotToRDF
are wired in so far. Add a source by giving it its own `rule <name>:` below
and listing it in `rule all`'s inputs.

plantTFDB's, stringdb's and uniprot's raw input is fetched on demand
(download_* rules) and deleted again as soon as it's converted (temp()
outputs) -- unlike interpro, which still expects data already placed under
data_dir by hand. Run with `snakemake --keep-going` so a species one of
these sources has no data for doesn't block the others.

Most rules are light enough that profiles/slurm/config.yaml's default
resources cover them; give a heavier future rule its own `resources:` block
(e.g. a bigger mem_mb_per_cpu or a different slurm_partition) rather than
raising that shared default.

Every source's Turtle output is validated with rapper before rule all ever
sees it, via the generic validate_turtle rule below -- a rule that produces
Turtle names its own output *.raw.ttl (wrapped in temp(), since
validate_turtle is its only consumer) instead of *.ttl, and validate_turtle
turns it into the final *.ttl by the naming convention alone, no per-rule
wiring needed. Catches malformed syntax from a parser bug (escaping,
encoding, ...) in one pass over the whole finished file, the same mechanism
ensembl_rdf's bin/rapper_postprocess.sh uses. This makes `rapper` (raptor2)
a hard dependency of every rule's output; validate_turtle loads it itself
(same defensive module-loading block as ensembl_rdf's bin/convert_slurm.sh
and bin/rapper_postprocess.sh use, overridable the same way via MODULES/
SITE_INIT) rather than relying on the submitting shell's own `module load`
having propagated to the job -- confirmed on the cluster not to, since a
non-interactive batch shell may never source the site's module init script.
"""
import os
from snakemake.exceptions import WorkflowError
from riceKB.utils import load_agrold_taxon_ids

configfile: "config.yaml"

# riceKB/*.py scripts do `from riceKB.globalVars import *`, which needs the
# repo root (not riceKB/ itself, the directory `python3 riceKB/x.py` puts on
# sys.path[0] by default) on PYTHONPATH. Don't rely on the submitting
# shell's own PYTHONPATH already having it (it won't, e.g. on a cluster
# compute node) -- set it explicitly for every rule's shell command instead.
shell.prefix(f'export PYTHONPATH="{workflow.basedir}:$PYTHONPATH"; ')

OUTPUT_DIR = config["output_dir"]

# {production_name: taxon_id} for every AgroLD species -- the single source
# of truth for which species exist and what each one's NCBI taxon is, so
# nothing that keys on either one (this Snakefile, bin/prepare_virtuoso_import.py)
# disagrees with the others. See riceKB/utils.py:load_agrold_taxon_ids.
try:
    AGROLD_TAXON_ID = load_agrold_taxon_ids(
        config["ensembl_rdf"]["species_file"], config["ensembl_rdf"]["species_table"])
except RuntimeError as e:
    raise WorkflowError(str(e))


# Like stringdb below, plantTFDB/PlantRegMap's raw input is downloaded on
# demand by the download_planttfdb rule, so every AgroLD species is a
# candidate; a species neither site covers simply fails its own
# download_planttfdb job (see bin/download_planttfdb.py) -- run with
# `snakemake --keep-going` so one such failure doesn't block the others.
PLANTTFDB_SPECIES = sorted(AGROLD_TAXON_ID)


# Unlike planttfdb (no fetch script wired in yet), stringdb's raw input is
# downloaded on demand by the download_stringdb rule below, so there is no
# pre-existing data_dir to scan: every AgroLD species is a candidate. A
# species STRING has no data for simply fails its own download_stringdb job
# (see bin/download_stringdb.py) -- run with `snakemake --keep-going` so one
# such failure doesn't block the others.
STRINGDB_SPECIES = sorted(AGROLD_TAXON_ID)


# UniProt publishes ONE flat-file per taxonomic division (e.g. "plants"),
# covering every species in it -- unlike planttfdb/stringdb above, there is
# no per-species variant to fetch, so this is keyed by taxon id rather than
# species name. Letting riceKB/uniprotToRDF.py filter its one taxon out of
# that shared file directly, once per taxon, meant re-decompressing and
# re-scanning the whole file (up to ~14GB compressed for trembl) once per
# taxon -- cheap for a handful of species, but that's what pushed some
# trembl taxa past even a generous SLURM runtime. split_uniprot below scans
# it exactly once instead, fanning every record out to whichever of
# UNIPROT_TAXA it belongs to; riceKB/uniprotToRDF.py itself didn't need to
# change, since it already reads plain text or gzip input the same way.
UNIPROT_BANKS = ["sprot", "trembl"]
UNIPROT_TAXA = sorted(set(AGROLD_TAXON_ID.values()))


rule all:
    input:
        f"{OUTPUT_DIR}/interpro.ttl",
        f"{OUTPUT_DIR}/oryzabase_genes.ttl",
        f"{OUTPUT_DIR}/reference_oryzabase.ttl",
        f"{OUTPUT_DIR}/qtaro_genes.ttl",
        f"{OUTPUT_DIR}/qtaro_qtl.ttl",
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/stringdb_{{species}}.ttl", species=STRINGDB_SPECIES),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.plants.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.keyword.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.pubmed.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/mappings_{{taxon}}_{{bank}}_gene_prot.tsv", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),


rule uniprot_all:
    """Just the uniprot outputs (both banks, every taxon) -- a convenience
    target for running/testing this one source without the rest of the
    pipeline's other data sources and their own prerequisites (e.g. a
    manually-placed data/oryzabase export that `all` also depends on).
    """
    input:
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.plants.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.keyword.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.pubmed.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/mappings_{{taxon}}_{{bank}}_gene_prot.tsv", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),


rule validate_turtle:
    """Validate and normalize a pipeline's raw Turtle output with rapper --
    the same mechanism ensembl_rdf's bin/rapper_postprocess.sh uses, applied
    here as one generic rule instead of a per-source step. {name} is a free
    wildcard (can contain dots, e.g. "uniprot_4530_trembl.plants") matched
    against whatever *.raw.ttl a producing rule's output already names.

    rapper's own exit code is 1 on a genuine syntax error (verified against
    a deliberately malformed file), so Snakemake's normal "remove a failed
    job's declared outputs" behaviour already keeps a bad {output} from
    lingering. {output}.rapper.log is kept only when something failed, for
    debugging -- on success its content is already visible in Snakemake's
    own job log via the cat below, so it's removed rather than left as
    permanent clutter next to every one of the pipeline's Turtle files.
    """
    input:
        f"{OUTPUT_DIR}/{{name}}.raw.ttl",
    output:
        f"{OUTPUT_DIR}/{{name}}.ttl",
    wildcard_constraints:
        # Without this, {name} is free to end in ".raw" itself, so Snakemake
        # considers applying this very rule to its own output recursively
        # (name="X.raw" producing "X.raw.ttl" as an input to name="X.raw.raw"
        # producing "X.raw.raw.ttl", ...) -- confirmed: Snakemake rejects the
        # DAG with a PeriodicWildcardError/AmbiguousRuleException without
        # this constraint.
        name=r".+(?<!\.raw)",
    shell:
        # Same defensive module-loading block as ensembl_rdf's
        # bin/convert_slurm.sh/bin/rapper_postprocess.sh (see module
        # docstring above for why it's needed): `module` is a shell
        # function a non-login batch shell may not have sourced, so this
        # loads the site's module init script itself before calling
        # `module load`, rather than assuming the submitting shell's own
        # `module load` reached this job. MODULES/SITE_INIT are
        # overridable from the environment exactly as in those scripts.
        "MODULES=${{MODULES-\"bioinfo-itrop raptor2/2.0.16\"}}; "
        "if [ -n \"$MODULES\" ]; then "
        "for f in \"${{SITE_INIT:-}}\" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh "
        "\"${{MODULESHOME:-}}/init/bash\" /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do "
        "[ -n \"$f\" ] && [ -f \"$f\" ] && . \"$f\" && break; done; "
        "if type module >/dev/null 2>&1; then "
        "for m in $MODULES; do module load \"$m\" || echo \"Warning: could not load module $m\" >&2; done; "
        "else echo \"Warning: no module command; skipping $MODULES\" >&2; fi; "
        "fi; "
        "rapper -i turtle -o turtle {input} > {output} 2> {output}.rapper.log; "
        "status=$?; cat {output}.rapper.log >&2; "
        "if [ $status -eq 0 ]; then rm -f {output}.rapper.log; fi; "
        "exit $status"


rule interpro:
    """InterPro domain hierarchy, labels and GO mapping."""
    input:
        tree_file=config["interpro"]["tree_file"],
        entry_list=config["interpro"]["entry_list"],
        go_mapping=config["interpro"]["go_mapping"],
    output:
        temp(f"{OUTPUT_DIR}/interpro.raw.ttl"),
    shell:
        "python3 riceKB/interpro.py {input.tree_file} {input.entry_list} "
        "{input.go_mapping} {output}"


rule oryzabase:
    """OryzaBase gene list -- names, symbols, RAP-DB/MSU/Gramene cross-refs,
    GO/TO/PO classifications -- plus the side mapping files other
    AgroLD_ETL parsers (see riceKB/rapdbModel.py) join against.

    OryzaBase's gene list is exported on demand via a web form
    (https://shigen.nig.ac.jp/rice/oryzabase/gene/download), with a
    timestamp baked into the downloaded filename -- no stable URL to fetch
    automatically, so like interpro above, it's placed by hand under
    data_dir (rename or symlink the download to the path config.yaml
    expects).
    """
    input:
        config["oryzabase"]["gene_list"],
    output:
        genes=temp(f"{OUTPUT_DIR}/oryzabase_genes.raw.ttl"),
        pub_dict=f"{OUTPUT_DIR}/pub_dictionnary.txt",
        pub_dict_rapdb=f"{OUTPUT_DIR}/pub_dictionnary_rapdb.txt",
        pub_dict_msu=f"{OUTPUT_DIR}/pub_dictionnary_msu.txt",
        pub_dict_hub=f"{OUTPUT_DIR}/pub_genome_hub.txt",
        dict_po=f"{OUTPUT_DIR}/dict_PO_hub.txt",
        dict_to=f"{OUTPUT_DIR}/dict_TO_hub.txt",
    shell:
        "python3 riceKB/oryzaBaseParser.py {input} {output.genes} "
        "--pub-dict {output.pub_dict} --pub-dict-rapdb {output.pub_dict_rapdb} "
        "--pub-dict-msu {output.pub_dict_msu} --pub-dict-hub {output.pub_dict_hub} "
        "--dict-po {output.dict_po} --dict-to {output.dict_to}"


rule oryzabase_reference:
    """Gene<->PubMed links from OryzaBase's separate "Reference" export
    (literature, not genes -- distinct download from oryzabase above),
    joined against the oryzabase rule's pub_dict output by gene
    symbol/synonym. See riceKB/oryzabaseReference.py.
    """
    input:
        reference=config["oryzabase"]["reference"],
        pub_dict=f"{OUTPUT_DIR}/pub_dictionnary.txt",
    output:
        temp(f"{OUTPUT_DIR}/reference_oryzabase.raw.ttl"),
    shell:
        "python3 riceKB/oryzabaseReference.py {input.reference} {input.pub_dict} {output}"


rule qtaro_genes:
    """QTARO's gene export. A handful of rows give an MSU locus id
    (LOC_Os...) instead of a RAP/IRGSP one; --msu-rap-mapping converts
    those via RAP-DB's RAP-MSU mapping (config["rapdb_msu"]) instead of
    dropping them, see riceKB/utils.py:load_msu_rap_mapping.

    Tracked in git under data/qtaro/ (small enough, see config.yaml), so
    unlike oryzabase above there's nothing to place by hand here.
    """
    input:
        gene_export=config["qtaro"]["gene_export"],
        mapping=config["rapdb_msu"]["mapping_file"],
    output:
        temp(f"{OUTPUT_DIR}/qtaro_genes.raw.ttl"),
    shell:
        "python3 riceKB/qtaroModel.py genes {input.gene_export} {output} "
        "--msu-rap-mapping {input.mapping}"


rule qtaro_qtl:
    """QTARO's QTL export (no MSU ids to convert here, unlike qtaro_genes
    above -- see config["qtaro"]["qtl_export"])."""
    input:
        config["qtaro"]["qtl_export"],
    output:
        temp(f"{OUTPUT_DIR}/qtaro_qtl.raw.ttl"),
    shell:
        "python3 riceKB/qtaroModel.py qtl {input} {output}"


rule download_planttfdb:
    """Fetch one species' PlantTFDB TF list + PlantRegMap merged regulation
    network on demand, via bin/download_planttfdb.py. Declared temp() so
    Snakemake deletes them as soon as the planttfdb rule below has consumed
    them -- same reasoning as download_stringdb below: no permanent local
    copy once converted.

    Unlike download_stringdb, output filenames here don't embed a taxon id
    (PlantTFDB/PlantRegMap files are always named the same within a species'
    subdirectory), so {species} is the only wildcard needed and no
    taxon-matching trick is required.
    """
    output:
        tf_list=temp(os.path.join(
            config["planttfdb"]["data_dir"], "{species}", config["planttfdb"]["tf_list"])),
        regulation_file=temp(os.path.join(
            config["planttfdb"]["data_dir"], "{species}", config["planttfdb"]["regulation_file"])),
    shell:
        "python3 bin/download_planttfdb.py -s {wildcards.species} "
        "-o {config[planttfdb][data_dir]} "
        "--tf-list-filename {config[planttfdb][tf_list]} "
        "--regulation-filename {config[planttfdb][regulation_file]} "
        "-f {config[ensembl_rdf][species_file]}"


rule planttfdb:
    """PlantTFDB gene families and the regulatory relations between them, for
    one species. The taxon id comes from the shared ensembl_rdf species
    list (load_agrold_taxon_ids above), not from the input file, which
    carries no species information of its own.

    One parser invocation writes both outputs; declaring them together as
    Snakemake output means a resume treats them as produced atomically (one
    missing file reruns the whole rule, not a partial one).
    """
    input:
        tf_list=lambda wc: os.path.join(
            config["planttfdb"]["data_dir"], wc.species, config["planttfdb"]["tf_list"]),
        regulation_file=lambda wc: os.path.join(
            config["planttfdb"]["data_dir"], wc.species, config["planttfdb"]["regulation_file"]),
    output:
        tf=temp(f"{OUTPUT_DIR}/planttfdb_{{species}}.raw.ttl"),
        regulation=temp(f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.raw.ttl"),
    params:
        taxon_id=lambda wc: AGROLD_TAXON_ID[wc.species],
    shell:
        "python3 riceKB/plantTFDB.py {input.tf_list} {output.tf} "
        "{input.regulation_file} {output.regulation} --taxon-id {params.taxon_id}"


rule download_stringdb:
    """Fetch one species' STRING aliases + links files on demand, via
    bin/download_stringdb.py. Declared temp() so Snakemake deletes them as
    soon as the stringdb rule below has consumed them: STRING's raw,
    unfiltered links file is several GB per species on disk, far more than
    the filtered RDF it gets converted to, so there is no reason to keep a
    permanent local copy per species once converted.

    {taxon} is a free wildcard here (not looked up from AGROLD_TAXON_ID
    directly, since output paths can't be computed by a function in
    Snakemake) -- it gets bound by matching this rule's output pattern
    against the concrete path the stringdb rule's input asks for, which
    already embeds the right taxon id.
    """
    output:
        aliases_file=temp(os.path.join(
            config["stringdb"]["data_dir"], "{species}",
            "{taxon}." + config["stringdb"]["aliases_suffix"])),
        links_file=temp(os.path.join(
            config["stringdb"]["data_dir"], "{species}",
            "{taxon}." + config["stringdb"]["links_suffix"])),
    wildcard_constraints:
        taxon=r"\d+",
    shell:
        "python3 bin/download_stringdb.py -s {wildcards.species} "
        "-o {config[stringdb][data_dir]} "
        "-f {config[ensembl_rdf][species_file]} "
        "-t {config[ensembl_rdf][species_table]}"


rule stringdb:
    """STRING protein-protein association network for one species, with
    STRING's own protein ids resolved to native gene locus ids via the
    aliases file (see riceKB/stringDB.py). Output is large -- for
    Arabidopsis thaliana alone the links file's ~15M rows produce a ~5GB
    Turtle file -- so this rule may need its own heavier `resources:` for
    bigger genomes, per the module docstring above.
    """
    input:
        aliases_file=lambda wc: os.path.join(
            config["stringdb"]["data_dir"], wc.species,
            f"{AGROLD_TAXON_ID[wc.species]}.{config['stringdb']['aliases_suffix']}"),
        links_file=lambda wc: os.path.join(
            config["stringdb"]["data_dir"], wc.species,
            f"{AGROLD_TAXON_ID[wc.species]}.{config['stringdb']['links_suffix']}"),
    output:
        temp(f"{OUTPUT_DIR}/stringdb_{{species}}.raw.ttl"),
    params:
        taxon_id=lambda wc: AGROLD_TAXON_ID[wc.species],
    shell:
        "python3 riceKB/stringDB.py {input.aliases_file} {input.links_file} "
        "{output} --taxon-id {params.taxon_id}"


rule download_uniprot:
    """Fetch one shared UniProt taxonomic division flat-file (one per bank:
    sprot or trembl), still gzip-compressed, via bin/download_uniprot.py --
    not per species, unlike stringdb/planttfdb above: UniProt publishes ONE
    file per division, covering every species in it, and
    riceKB/uniprotToRDF.py filters the species it wants out of this shared
    file itself (see the uniprot rule below). Kept gzip-compressed on disk
    (not decompressed here) because uniprot_trembl_plants.dat.gz is ~14GB
    compressed -- decompressing it to a second on-disk copy isn't viable on
    a quota-limited cluster; riceKB/uniprotToRDF.py decompresses on the fly
    as it reads instead. Declared temp() so Snakemake deletes it once every
    uniprot job for this bank, across every taxon, has consumed it.
    """
    output:
        temp(os.path.join(config["uniprot"]["data_dir"], "uniprot_{bank}_plants.dat.gz")),
    wildcard_constraints:
        bank=r"sprot|trembl",
    shell:
        "python3 bin/download_uniprot.py --bank {wildcards.bank} --division plants "
        "-o {config[uniprot][data_dir]}"


rule split_uniprot:
    """Split one bank's shared plants division dump into one plain-text
    file per taxon in UNIPROT_TAXA, in a single pass (bin/
    split_uniprot_by_taxon.py), instead of letting the uniprot rule below
    re-decompress and re-scan the whole (up to ~14GB compressed, for
    trembl) file once per taxon -- see UNIPROT_TAXA's comment above for why
    that mattered. Always writes a file for every requested taxon, even an
    empty one for a taxon this bank has no entry for, so a taxon with
    nothing to convert still fails the uniprot job below exactly as it did
    reading directly from the undivided dump (see that rule's docstring).

    One job per bank, not per taxon: declared temp() per output taxon file,
    so Snakemake deletes each one once its own uniprot job (not every
    taxon's -- each taxon has its own file here) has consumed it.
    """
    input:
        os.path.join(config["uniprot"]["data_dir"], "uniprot_{bank}_plants.dat.gz"),
    output:
        [temp(f"{OUTPUT_DIR}/uniprot_split/{{bank}}/{taxon}.dat") for taxon in UNIPROT_TAXA],
    wildcard_constraints:
        bank=r"sprot|trembl",
    params:
        output_dir=f"{OUTPUT_DIR}/uniprot_split/{{bank}}",
        taxa=" ".join(UNIPROT_TAXA),
    resources:
        # The profile's 30-minute default was sized for the small per-rule
        # jobs elsewhere in this pipeline; this rule's single pass over the
        # full trembl dump (~14GB compressed) is the one job that used to
        # time out when done per-taxon (see UNIPROT_TAXA's comment above),
        # so give it room even though it now runs just once per bank.
        runtime=360,
    shell:
        "python3 bin/split_uniprot_by_taxon.py {input} {params.output_dir} {params.taxa}"


rule uniprot:
    """UniProt protein records -- plus the keywords and PubMed citations
    they reference -- for one taxon and one bank (sprot or trembl), out of
    that taxon's own split (split_uniprot above) of the plants division
    dump.

    {taxon} is a free wildcard here (not looked up via AGROLD_TAXON_ID,
    since output paths can't be computed by a function in Snakemake); it's
    rule all's expand() over AGROLD_TAXON_ID's values (UNIPROT_TAXA above)
    that actually drives which taxa get requested. A taxon with no entry
    for this bank in the plants division fails this job (see
    riceKB/uniprotToRDF.py's exit(1) on zero matches) -- run with
    `snakemake --keep-going` so that doesn't block the other taxa/banks.

    riceKB/uniprotToRDF.py builds its own output filenames internally from
    bank/taxon (uniprot_{taxon}_{bank}.plants.ttl etc. -- its CLI doesn't
    take explicit output paths, deliberately left unchanged), so unlike the
    other sources this rule can't just rename its own `output:` to
    *.raw.ttl and let the script write there directly; it renames the
    script's fixed-name files to *.raw.ttl itself afterwards instead, so
    validate_turtle still picks them up the same way as everywhere else.
    """
    input:
        f"{OUTPUT_DIR}/uniprot_split/{{bank}}/{{taxon}}.dat",
    output:
        protein=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.plants.raw.ttl"),
        keyword=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.keyword.raw.ttl"),
        pubmed=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.pubmed.raw.ttl"),
        gene_mapping=f"{OUTPUT_DIR}/mappings_{{taxon}}_{{bank}}_gene_prot.tsv",
    wildcard_constraints:
        taxon=r"\d+",
        bank=r"sprot|trembl",
    params:
        output_dir=OUTPUT_DIR,
    shell:
        "python3 riceKB/uniprotToRDF.py {wildcards.bank} {wildcards.taxon} "
        "{params.output_dir} {input} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.plants.ttl {output.protein} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.keyword.ttl {output.keyword} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.pubmed.ttl {output.pubmed}"
