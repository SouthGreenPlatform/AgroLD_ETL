
'''
This file contains all global variables i.e. URI patterns, Ontological terms and predicates to be used by the various parsers
 
'''
import re
'''
TODO:
    1) Add variables for common predicates used in the Rdf converters
    2) Cleanup cleanup the script and add comments
'''
global taxon_ids, db_obj_type, base, pr, rdf, rdf_ns, rdfs_ns, skos, skos_ns, dc_ns, dc_uri, dcterms_ns, dcterms_uri, faldo, faldo_ns,\
    ncbi_tax_ns, ncbi_tax_uri, owl, owl_ns, xsd, xsd_ns, base_uri, base_ns, base_vocab_uri, base_vocab_ns,\
    sio_uri, sio_ns, obo_uri, obo_ns, ncbi_tax_uri, ncbi_tax_ns, uniprot, up_ns, interpro_uri, interpro_ns, chromosome_ns, chromosome_uri,\
    ensembl_gene_uri, ensembl_gene_ns, ensembl_transcript_uri, ensembl_transcript_ns, ensembl_protein_uri, ensembl_protein_ns,\
    ensembl_ns, ensembl_plant, gramene_gene, gr_g_ns, gramene_qtl, gr_qtl_ns, sio_term, go_aspects, ont_aspects,\
    gene_term, protein_term, tigr_uri, tigr_ns, identifiers_uri, rapdb_gene_uri, rapdb_gene_ns, rapdb_mrna_ns, rapdb_mrna_uri,\
    plant_trait_term,orygene_uri, orygene_ns, goa_uri, goa_ns, gr_assoc, gr_assoc_ns, tair_l_uri, tair_l_ns,\
    met_pw_sio_term, ec_code_uri, ec_code_ns,reaction_uri, reaction_ns, pathway_uri, pathway_ns, otl_uri, otl_ns,\
    plant_dev_term, plant_anatomy_term,germplasm_term, co_uri, co_ns, swo_uri, swo_ns, biocyc_pw_term, biocyc_react_term,\
    string_ns, string_uri, rap_pattern, gramene_pattern, prot_pattern, tigr_pattern, tair_pattern, ont_pattern,\
    string_pattern, sorghum_pattern, alt_sorghum_match, arabidopsis_pattern, maize_pattern, alt_maize_match,\
    pubmed_pattern, ncbi_pattern, interpro_pattern, pfam_pattern, ena_embl_uri, ena_embl_ns,\
    base_resource_uri, base_resource_ns, up_core_uri, up_core_ns, bibo_uri, bibo_ns, prism_uri, prism_ns

# Taxon - 'NCBI taxon IDs' : 'Taxon name' 
taxon_ids = {
        '4615':'Ananas comosus',
        '13333':'Amborella trichopoda',
        '3702': 'Arabidopsis thaliana',
        '3818':'Arachis hypogaea',
        '3708': 'Brassica napus',
        '3712':'Brassica oleracea',
        '1753':'Brassica_rapa',
        '3821':'Cajanus_cajan',
        '4543': 'Cenchrus americanus',
        '3827':'Cicer arietinum',
        '85681':'Citrus clementina',
        '49390': 'Coffea canephora str. DH200-94',
        '3847' :'Glycine max',
        '29730' :'Gossypium raimondii',
        '4513' : 'Hordeum vulgare',
        '3983': 'Manihot esculenta',
        '3880' : 'Medicago truncatula',
        '214687': 'Musa acuminata subsp. malaccensis',
        '49451' : 'Nicotiana attenuata',
        '4543' : 'Cenchrus americanus',
        '3885' : 'Phaseolus vulgaris',
        '65489': 'Oryza barthii',
        '4533' : 'Oryza brachyantha',
        '4538' : 'Oryza glaberrima',
        '40148' : 'Oryza_glumaepatula',
        '4528': 'Oryza longistaminata',
        '40149': 'Oryza meridionalis',
        '4536': 'Oryza_nivara',
        '4537': 'Oryza_punctata',
        '4529': 'Oryza rufipogon',
        '4530' : 'Oryza sativa',
        '39946': 'Oryza indica',
        '39947' : 'Oryza sativa',
        '3694' : 'Populus trichocarpa',
        '62335' : 'Saccharum spontaneum',
        '4555'  : 'Setaria italica',
        '4558' : 'Sorghum bicolor',
        '3641' : 'Theobroma cacao',
        '4565' : 'Triticum aestivum',
        '85692' : 'Triticum dicoccoides',
        '4571'  :  'Triticum turgidum',
        '4572' : 'Triticum urartu',
        '29760' : 'Vitis vinifera',
        '4577' : 'Zea mays',
         }


# RegEX patterns
rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
gramene_pattern = re.compile(r'^GR\:\d{7}$')
prot_pattern = re.compile(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')
tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.*\d*$')
tair_pattern = re.compile(r'^AT[1-5]G\d{5}$')
ont_pattern = re.compile(r'^\w+\:\d{7}$')
string_pattern = re.compile(r'^([A-N,R-Z][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9])|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])|([0-9][A-Za-z0-9]{3})$')
sorghum_pattern = re.compile(r'^SB\d{2}G\d{6}\.\d$')
alt_sorghum_match = re.compile(r'^SB\d{4}S\d{6}\.\d$')
arabidopsis_pattern = re.compile(r'^AT[1-5]G\d{5}$')
maize_pattern = re.compile(r'^GRMZM\d{1}G\d{6}')
alt_maize_match = re.compile(r'^\w+\d{6}\.\d{1}\_\w+\d{3}')
pubmed_pattern = re.compile(r'^\d+$')
ncbi_pattern = re.compile(r'^[A-Z]{2}\d{6}$')
interpro_pattern = re.compile(r'^IPR[0-9]{6}$')
pfam_pattern = 	re.compile(r'^PF\d{5}$')

# Resolvable URIs
db_obj_type = {
              'protein' : 'http://www.identifiers.org/uniprot/', # Note: assumes protein accessions are from UniProt
              'gene' : 'http://www.identifiers.org/gramene.gene/', # Note: assumes gene accessions are Gramene internal IDs e.g. GR:xxxxx
              'QTL' : 'http://www.identifiers.org/gramene.qtl/' # Note: assumes QTL accessions are Gramene QTL IDs e.g. AQED049
              }

# Prefixes
base = '@base'
base_uri = 'http://purl.agrold.org/'
base_ns = 'agrold:'

pr = '@prefix'

rdf = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
rdf_ns = 'rdf:'

rdfs = 'http://www.w3.org/2000/01/rdf-schema#'
rdfs_ns = 'rdfs:'

owl = 'http://www.w3.org/2002/07/owl#'
owl_ns = 'owl:'

xsd = 'http://www.w3.org/2001/XMLSchema#'
xsd_ns = 'xsd:'

skos = 'http://www.w3.org/2004/02/skos/core#'
skos_ns = 'skos:'

owl_uri = 'http://www.w3.org/2002/07/owl#'
owl_ns = 'owl:'

dc_uri = 'http://purl.org/dc/elements/1.1/'
dc_ns = 'dc:'

dcterms_uri = 'http://purl.org/dc/terms/'
dcterms_ns = 'dcterms:'

doi_uri = 'http://dx.doi.org/'
doi_ns = 'doi:'

bibo_uri = 'http://purl.org/ontology/bibo/'
bibo_ns = 'bibo:'

prism_uri = 'http://prismstandard.org/namespaces/1.2/basic/'
prism_ns = 'prism:'
#Internal URI/namespaces
#base_uri = 'http://purl.agrold.org/resource/'
#base_ns = 'agrold:'

base_vocab_uri = 'http://purl.agrold.org/vocabulary/'
base_vocab_ns = 'agrold_vocabulary:'

base_resource_uri = 'http://purl.agrold.org/resource/'
base_resource_ns = 'resource:'
# Datasource specific URIs
ncbi_gene_uri = 'http://identifiers.org/ncbigene/'
ncbi_gene_ns = 'ncbigene:'

ncbi_tax_uri = 'http://identifiers.org/taxonomy/'
ncbi_tax_ns = 'taxon:'

#uniprot = 'http://www.identifiers.org/uniprot/'
uniprot = 'http://purl.uniprot.org/uniprot/'
up_ns = 'uniprot:'

faldo = 'http://biohackathon.org/resource/faldo#'
faldo_ns = 'faldo:'

gramene_gene =  'http://www.identifiers.org/gramene.gene/'
gr_g_ns = 'gramene_gene:'

resource = 'http://purl.agrold.org/resource/'
res_ns = 'agrold_resource:'

gramene_qtl = 'http://www.identifiers.org/gramene.qtl/'
gr_qtl_ns = 'gramene_qtl:' 

msu_uri = 'http://purl.agrold.org/resource/msu.locus/'
msu_ns = 'msu:'

qtaro_qtl = 'http://purl.agrold.org/resource/qtaro.qtl/'
qtaro_qtl_ns = 'qtaro_qtl:'


ensembl_plant = 'http://identifiers.org/ensembl.plant/' # http://rdf.ebi.ac.uk/resource/ensembl/
ensembl_ns = 'ensembl:'

ensembl_gene_uri = 'http://rdf.ebi.ac.uk/resource/ensembl/'
ensembl_gene_ns = 'ensembl_gene:'


ensembl_transcript_uri = 'http://rdf.ebi.ac.uk/resource/ensembl.transcript/'
ensembl_transcript_ns = 'ensembl_transcript:'


ensembl_protein_uri = 'http://rdf.ebi.ac.uk/resource/ensembl.protein/'
ensembl_protein_ns = 'ensembl_protein:'


tair_l_uri = 'http://identifiers.org/tair.locus/'
tair_l_ns = 'tairlocus:'

tigr_uri = 'http://purl.agrold.org/resource/msu/'
tigr_ns = 'tigr:'

tigr_g_uri = 'http://identifiers.org/ricegap/'
tigr_g_ns = 'tigr_gene:'

rapdb_gene_uri = 'http://identifiers.org/rapdb.locus/'
rapdb_gene_ns = 'rapdb_gene:'

rapdb_mrna_uri = 'http://identifiers.org/rapdb.mrna/'
rapdb_mrna_ns = 'rapdb_mrna:'


tenor_uri = 'http://tenor.dna.affrc.go.jp/EPV/'
tenor_ns = 'tenor:'

interpro_uri = 'http://purl.uniprot.org/interpro/'
interpro_ns = 'interpro:'

uniprot_uri = 'http://identifiers.org/uniprot/'
uniprot_ns = 'uniprot:'

up_base_uri = 'http://purl.uniprot.org/'
up_base_ns = 'up_base_ns:'

up_core_uri = 'http://purl.uniprot.org/core'
up_core_ns = 'up_core_ns:'

orygene_ns = 'http://identifiers.org/oryzabase.gene/'
orygene_uri = 'oryzabase:'

oryzabase_uri = 'http://identifiers.org/oryzabase.gene/'
oryzabase_ns = 'oryzabase:'

ec_code_uri = 'http://identifiers.org/ec-code/'
ec_code_ns = 'ec:'

reaction_uri = 'http://purl.agrold.org/resource/biocyc.reaction/'
reaction_ns = 'reaction:'

pathway_uri = 'http://purl.agrold.org/resource/biocyc.pathway/'
pathway_ns = 'pathway:'

#BioCyc
swo_uri = 'http://edamontology.org/'
swo_ns = 'swo:'

#Pubmed
pubmed_uri = 'http://identifiers.org/pubmed/'
pubmed_ns = 'pubmed:'

pubmed_rdf_uri = 'http://rdf.ncbi.nlm.nih.gov/pubmed/'
pubmed_rdf_ns = 'pubmed_rdf:'
# AraCyc
#aracyc_uri = 'http://purl.agrold.org/resource/aracyc.pathway/'
#aracyc_ns = 'aracyc_pathway:'
#aracyc_gene_uri = 'http://purl.agrold.org/resource/aracyc.gene/'
#aracyc_gene_ns = 'aracyc_gene:'
#aracyc_prot_uri = 'http://purl.agrold.org/resource/aracyc.protein/'
#aracyc_prot_ns = 'aracyc_protein:'

#RiceCyc
#ricecyc_uri = 'http://purl.agrold.org/resource/ricecyc.pathway/'
#ricecyc_ns = 'ricecyc_pathway:'

# SouthGreen
#----------------
# OTL
otl_uri = 'http://identifiers.org/otl/'
otl_ns = 'otl:'

# identifiers
identifiers_uri = 'http://identifiers.org/'


#Orygenesdb
mirbase_uri = 'http://www.identifiers.org/mirbase/'
mirbase_ns = 'mirbase:'

mirbase2_uri ='http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc='
mirbase2_ns = 'mirbase:'

mirbase_mature_uri = 'http://www.identifiers.org/mirbase.mature/'
mirbase_mature_ns = 'mature:'

protein_uri = 'http://purl.agrold.org/resource/orygenesdb.protein/'
protein_ns = 'protein:'

flanking_uri = 'http://purl.agrold.org/resource/orygenesdb.flankingregion/'
flanking_ns = 'flanking:'

mRNA_uri = 'http://purl.agrold.org/resource/orygenesdb.mrna/'
mRNA_ns = 'mrna:'

cDNA_uri = 'http://purl.agrold.org/resource/orygenesdb.cdna/'
cDNA_ns = 'cdna:'

otl_public_plante_uri = 'http://www.identifiers.org/otl/'
otl_public_plante_ns = 'otl_plante:'

OrygenesDB_uri = 'http://purl.agrold.org/resource/orygenesdb/'
OrygenesDB_ns = 'orygenesdb:'

RiceNetDB_gene_uri = 'https://identifiers.org/ricenetdb.gene/'
RiceNetDB_gene_ns = 'ricenetdb_gene:'


RiceNetDB_protein_uri = 'https://identifiers.org/ricenetdb.protein/'
RiceNetDB_protein_ns = 'ricenetdb_protein:'

chromosome_uri = 'http://purl.agrold.org/resource/chromosome/'
chromosome_ns = 'chromosome:'

ena_embl_uri = 'http://identifiers.org/ena.embl/'
ena_embl_ns = 'ena.embl:'

marker_uri = 'http://purl.agrold.org/resource/orygenesdb.marker/'
marker_ns = 'marker:'

kegg_uri='http://identifiers.org/kegg/'
kegg_ns='kegg:'


kegg_path_uri='http://identifiers.org/kegg.pathway/'
kegg_path_ns= 'kegg_path:'

metacyc_uri='https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object='
metacyc_ns='metacyc:'

osa_uri='http://www.genome.jp/dbget-bin/www_bget?'
osa_ns='osa:'

string_uri = 'http://identifiers.org/string:'
string_ns = 'string:'

dosa_uri='http://www.genome.jp/dbget-bin/www_bget?'
dosa_ns='kegg_gene:'
#TROPGENE

study_uri = 'http://purl.agrold.org/resource/tropgene.study/'
study_ns = 'study:'

population_uri = 'http://purl.agrold.org/resource/tropgene.population/'
population_ns = 'population:'

qtl_uri = 'http://purl.agrold.org/resource/tropgene.qtl/'
qtl_ns = 'qtl:'

edam_uri ='http://edamontology.org/'
edam_ns = 'edam:'

trait_uri = 'http://purl.agrold.org/resource/tropgene.trait/'
trait_ns = 'trait:'

mapfeature_uri = 'http://purl.agrold.org/resource/tropgene.mapfeature/'
mapfeature_ns = 'mapfeature:'



# Ontology terms and aspects
sio_term = 'SIO_000897' # SIO term - association
gene_term = 'SO_0000704' # SO term - gene
mrna_term = 'SO_0000234' # SO term - mrna
protein_term = 'SO_0000104' # SO term - protein
qtl_term = 'SO_0000771' # SO term - qtl
plant_trait_term = 'TO_0000387' # TO term - plant trait (root term)
plant_dev_term = 'PO_0009012' # PO term - plant_structure_development_stage (root term)
plant_anatomy_term = 'PO_0025131' # PO term - plant_anatomy (root term)
germplasm_term = 'CO_715:0000225' # CO term - Passport information entity
met_pw_sio_term = 'SIO_010532' # SIO term - metabolic pathway
biocyc_pw_term = 'data_1157' # Pathway ID (BioCyc)
biocyc_react_term = 'data_2106' # Reaction ID (BioCyc)

co_uri = 'http://www.cropontology.org/rdf/'
co_ns = 'co:'

sio_uri = 'http://semanticscience.org/resource/' # association URI
sio_ns = 'sio:'

obo_uri = 'http://purl.obolibrary.org/obo/'
obo_ns = 'obo:'

#goa_uri = 'http://purl.agrold.org/resource/go.association/'
goa_uri = 'http://identifiers.org/goa/' 
goa_ns = 'goa:'

gr_assoc = 'http://purl.agrold.org/resource/gramene.association/'
gr_assoc_ns = 'gramene_association:'



go_aspects = {           
           # GO aspects
           'P' : 'BFO_0000056', # participates_in -  Biological process
           'F' : 'BFO_0000085', # has_function -  Molecular Function
           'C' :  'BFO_0000082' # located_in - Cellular Component
           }
ont_aspects = {
           # GO aspects
           'P' : 'participates_in', # participates_in - BFO_0000056  Biological process
           'F' : 'has_function', # has_function - BFO_0000085  Molecular Function
           'C' :  'located_in', # located_in - BFO_0000082 Cellular Component           
           # PO aspects - temporary, suffixed with agrold: relative URI  
           'A' : 'expressed_in',
           'G' : 'expressed_at',
           # TO aspect - temporary
           'T' : 'has_trait',
           # EO aspects - temporary
#           'E' : 'has_condition'
           'E' : 'observed_in'
           }


#sniplay
sniplay_gene_uri = 'http://identifiers.org/ricegap/'
sniplay_gene_ns = 'ricegap:'

sniplay_gene_integenic_uri =  'http://purl.agrold.org/resource/sniplay.gene/'
sniplay_gene_integenic_ns = 'sniplay_gene:'

sniplay_individual_uri =  'http://purl.agrold.org/resource/sniplay.individual/'
sniplay_individual_ns = 'sniplay_individual:'

sniplay_allele_uri = 'http://purl.agrold.org/resource/sniplay.allele/'
sniplay_allele_ns = 'allele:'

sniplay_pos_uri =  'http://purl.agrold.org/resource/sniplay.pos/'
sniplay_pos_ns = 'sniplay_pos:'

sniplay_consequence_uri =  'http://purl.agrold.org/resource/sniplay.consequence/'
sniplay_consequence_ns = 'consequence:'


#Brapi urgi phenotype
observationUnitDbId_uri = 'http://purl.agrold.org/resource/phenotype.observationUnitDbId/'
observationUnitDbId_ns = 'observationUnitDbId:'

studyDbId_uri = 'http://purl.agrold.org/resource/phenotype.studyDbId/'
studyDbId_ns = 'studyDbId:'

treatment_uri = 'http://purl.agrold.org/resource/phenotype.treatment/'
treatment_ns = 'treatment:'

observationUnitLevelLabel_uri = 'http://purl.agrold.org/resource/phenotype.observationLabel/'
observationUnitLevelLabel_ns = 'observation_level_label:'

observationVariableDbId_uri = 'http://purl.agrold.org/resource/phenotype.observationData/'
observationVariableDbId_ns = 'observation_data:'


#Wine dataset file name :

repeatmasker_uri = 'http://purl.agrold.org/resource/winedb.repeatmasker/'
repeatmasker_ns = 'repeatmasker:'