'''
Updated on Dec 4, 2019

@author: larmande
'''
# TODO chercher dans les synonyms les pattern RAPDB et MSU pour creer une nouvelle relation
# TODO ajouter Prot:uri encodedBy Gene:uri
from Bio import SwissProt

from riceKB.globalVars import *
#from globalVars import *
from riceKB.utils import *
import pprint
import re
import os
import sys
from collections import defaultdict

def keyword2URI(keyword):
    cleanKey = re.sub('^\s+|\s+$','',keyword)
    cleanKey = re.sub('\s','_',keyword)
    keywordURI =  "<" + base_resource_uri + 'keyword/' + cleanKey + ">"
    return keywordURI, cleanKey

def keyword2Triples(cleanedKey, keywordBuffer):
    # cleanKey = re.sub('^\s+|\s+$', '', keyword)
    keyword = re.sub('_', ' ', cleanedKey)
    keywordBuffer += "<" + base_resource_uri + "keyword/" + cleanedKey + "> \n"
    keywordBuffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Keyword" + " ;\n"
    keywordBuffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (keyword) + " .\n"
    return keywordBuffer

def cleanUp(text, title, provenance=False):
    clean_text = text.replace('"', '')
    clean_text = clean_text.replace(title, '')
    clean_text = re.sub('^\s+', '', clean_text)
    if provenance:
        clean_text = re.sub('\{.+?\}', '', clean_text)
    return clean_text

def splitComments(comments):
    buffer = ''
    for annotation in comments:
        if 'FUNCTION:' in annotation:
            # (RO_0000085  has function)
            function = cleanUp(annotation,'FUNCTION:')
            buffer += "\t" + obo_ns + "RO_0000085" + "\t" + '"%s"' % (function) + " ;\n"
            # print('found FUNCTION:')
        elif 'SUBCELLULAR LOCATION:' in annotation:
            # RO_0001025  located in
            location = cleanUp(annotation,'SUBCELLULAR LOCATION:', True)
            if 'Note=' in location:
                location = location.split('Note=')[0]
            location_entries = re.split('\.|,|;',location)
            for location_entry in location_entries:
                location_entry = re.sub('^\s+|\s+$', '', location_entry)
                if location_entry:
                    buffer += "\t" + obo_ns + "RO_0001025" + "\t" + '"%s"' % (location_entry) + " ;\n"
            print('found SUBCELLULAR LOCATION:')
        elif 'TISSUE SPECIFICITY:' in annotation:
            #  up:Tissue_Specificity_Annotation RO_0002206  expressed in
            tissue = cleanUp(annotation,'TISSUE SPECIFICITY:')
            buffer += "\t" + obo_ns + "RO_0002206" + "\t" + '"%s"' % (tissue) + " ;\n"
            # Expressed in a small number of epidermal or subepidermal cells at the leaf axils, in axillary
            # meristems and the entire tiller buds.
            # Undetected in the shoot apical meristem. {ECO:0000269|PubMed:12687001}.
            print('found TISSUE SPECIFICITY:')
        elif 'DOMAIN:'  in annotation:
            #  RO_0002524:contains
            domain = cleanUp(annotation,'DOMAIN:')
            buffer += "\t" + obo_ns + "RO_0002524" + "\t" + '"%s"' % (domain) + " ;\n"
            # The C-terminal part of the protein is important for tillering.
            # Mutant moc1, in which the last 124 amino acids are missing, is mono culm. {ECO:0000305}.
            print('found DOMAIN:')
        elif 'SIMILARITY:' in annotation:
            #  Belongs to the GRAS family. {ECO:0000305}.
            print('found SIMILARITY:')
        elif 'SEQUENCE CAUTION:' in annotation:
            # Sequence=BAD35485.1; Type=Erroneous gene model prediction; Evidence={ECO:0000305};
            # Sequence=BAS98577.1; Type=Erroneous gene model prediction; Evidence={ECO:0000305};
            print('found SEQUENCE CAUTION:')
        elif 'CATALYTIC ACTIVITY:' in annotation:
#             ''' Reaction=cytidine(32)/guanosine(34) in tRNA + 2 S-adenosyl-L-methionine
# CC         = 2'-O-methylcytidine(32)/2'-O-methylguanosine(34) in tRNA + 2 H(+) +
# CC         2 S-adenosyl-L-homocysteine; Xref=Rhea:RHEA:42396, Rhea:RHEA-
# CC         COMP:10246, Rhea:RHEA-COMP:10247, ChEBI:CHEBI:15378,
# CC         ChEBI:CHEBI:57856, ChEBI:CHEBI:59789, ChEBI:CHEBI:74269,
# CC         ChEBI:CHEBI:74445, ChEBI:CHEBI:74495, ChEBI:CHEBI:82748; '''
            print('found CATALYTIC ACTIVITY:')
        elif 'COFACTOR:'  in annotation:
            #     Name=Zn(2+); Xref=ChEBI:CHEBI:29105;
            # CC         Evidence={ECO:0000269|PubMed:11076507, ECO:0000269|PubMed:12499545,
            # CC         ECO:0000269|PubMed:1336460, ECO:0000269|PubMed:1433293,
            # CC         ECO:0000269|PubMed:1909891, ECO:0000269|PubMed:19583303,
            # CC         ECO:0000269|PubMed:3151019, ECO:0000269|PubMed:3151020,
            # CC         ECO:0000269|PubMed:7761440, ECO:0000269|PubMed:7803386,
            # CC         ECO:0000269|PubMed:7901850, ECO:0000269|PubMed:8218160,
            # CC         ECO:0000269|PubMed:8262987, ECO:0000269|PubMed:8331673,
            # CC         ECO:0000269|PubMed:8399159, ECO:0000269|PubMed:8431430,
            # CC         ECO:0000269|PubMed:8451242, ECO:0000269|PubMed:8482389,
            # CC         ECO:0000269|PubMed:8639494, ECO:0000269|PubMed:8987974,
            # CC         ECO:0000269|PubMed:9398308, ECO:0000269|PubMed:9865942};
            # CC       Name=Co(2+); Xref=ChEBI:CHEBI:48828;
            # CC         Evidence={ECO:0000269|PubMed:19583303};
            # CC       Note=Zinc. Can also use cobalt(II) with lower efficiency, but not
            # CC       copper(II), nickel(II) and manganese(II).
            # CC       {ECO:0000269|PubMed:19583303};
            print('found COFACTOR:')
        elif 'ACTIVITY REGULATION:'  in annotation:
            # Phosphorylation leads to an increase in the catalytic activity.
            print('found ACTIVITY REGULATION:')
        elif 'BIOPHYSICOCHEMICAL PROPERTIES:'  in annotation:
            #            Kinetic parameters:
            # CC         KM=98 uM for ATP;
            # CC         KM=688 uM for pyridoxal;
            # CC         Vmax=1.604 mmol/min/mg enzyme;
            # CC       pH dependence:
            # CC         Optimum pH is 6.0. Active from pH 4.5 to 10.5.;
            print('found BIOPHYSICOCHEMICAL PROPERTIES:')
        elif 'SUBUNIT:' in annotation:
            # subPropertyOf Interaction(s)
            # Heterotrimer of alpha, beta and gamma subunits.
            # CC       {ECO:0000269|PubMed:12730181}.
            #  Self-associates. Interacts with BNIP3 and STEAP3. Interacts
            # CC       (via BH3 domain) with SPATA18 (via coiled-coil domains).
            # CC       {ECO:0000269|PubMed:10381623, ECO:0000269|PubMed:12606722,
            # CC       ECO:0000269|PubMed:21264228}.
            print('found SUBUNIT:')
        elif 'PATHWAY:'  in annotation:
            #  serveral PATHWAY: PATHWAY: PATHWAY: possible
            # Nitrogen metabolism; (S)-allantoin degradation.
            # CC       {ECO:0000255|HAMAP-Rule:MF_00616}.
            #  Each process is split up into 'super-pathway', 'pathway' and/or 'sub-pathway' associated with the protein.
            print('found PATHWAY:')
        elif 'ACTIVITY REGULATION:'  in annotation:
            print('found ACTIVITY REGULATION:')
        elif 'DEVELOPMENTAL STAGE:' in annotation:
            #  large block of text ..
            print('found DEVELOPMENTAL STAGE:')
        elif 'INDUCTION:' in annotation:
            #  By wounding, drought and salt stresses, benzothiadiazole
            # CC       (BTH), ethephon, hydrogen peroxide, abscisic acid (ABA) and
            # CC       incompatible and compatible races of rice blast fungus (M.grisea) and
            # CC       rice bacterial blight (X.oryzae). {ECO:0000269|PubMed:16766513}.
            print('found INDUCTION:')
        elif 'PTM:' in annotation:
            #  N-glycosylated and probably also O-glycosylated.
            print('found PTM:')
        elif 'RNA EDITING:' in annotation:
            # long text
            print('found RNA EDITING:')
        elif 'MASS SPECTROMETRY:' in annotation:
            # long text
            print('found MASS SPECTROMETRY:')
        elif 'POLYMORPHISM:' in annotation:
            #  long text
            print('found POLYMORPHISM:')
        elif 'DISEASE:' in annotation:
            # sseveral block  DISEASE:  DISEASE:
            #  Disease name (Disease abbreviation) [Link to OMIM]: Disease description.
            print('found  DISEASE:')
        elif 'DISRUPTION PHENOTYPE:' in annotation:
            #  text
            print('found DISRUPTION PHENOTYPE:')
        elif 'ALLERGEN:' in annotation:
            # Causes an allergic reaction in human.
            print('found ALLERGEN:')
        elif 'INTERACTION:' in annotation:
            # Q8K1M6; Q925I1: Atad3; NbExp=13; IntAct=EBI-2365792, EBI-772703;
            # Q8K1M6; Q5S006: Lrrk2; NbExp=5; IntAct=EBI-2365792, EBI-2693710;
            print('found INTERACTION:')
        elif 'TOXIC DOSE:' in annotation:
            # PD(50) is 0.36 mg/kg by injection in blowfly larvae, in
            # CC       PubMed:12885226.
            print('found TOXIC DOSE:')
        elif 'BIOTECHNOLOGY:' in annotation:
            # text
            print('found BIOTECHNOLOGY:')
        elif 'PHARMACEUTICAL:' in annotation:
            # text
            print('found PHARMACEUTICAL:')
        elif 'MISCELLANEOUS:' in annotation:
            # text
            print('found MISCELLANEOUS:')
        elif 'SEQUENCE CAUTION:' in annotation:
            # text
            print('found SEQUENCE CAUTION:')
        print(annotation)
    return buffer
def upToRDF(up_files, rdf_out_dir, additional_file):  # , output_file

    rdf_file = "uniprot.plants.ttl"
    rdf_keyword_file = "uniprot.keyword.ttl"
    output_file = os.path.join(rdf_out_dir, rdf_file)
    output_key_file = os.path.join(rdf_out_dir, rdf_keyword_file)
    output_writer = open(output_file, "w")
    keyword_writer = open(output_key_file, "w")
    keywordBuffer = ''
    rdf_buffer = ''
    prot_counter = 0
    lookup_list = set()
    keyword_list = set()
    pp = pprint.PrettyPrinter(indent=4)
    up_base_uri = "http://purl.uniprot.org/"
    #    up_base_ns = "uniprot_base:"
    if (additional_file):
        if (os.path.isfile(str(additional_file))):
            with open(additional_file, 'r') as infile:
                for prot in infile:
                    prot = re.sub('\s+', '', prot)
                    lookup_list.add(prot)

    print("************* Converting Uniprot data to RDF ***************\n")

    output_writer.write(str(getRDFHeaders()))
#    for upfile in up_files:
#        file_handle = open(upfile, "r")
    with open(up_files, 'r') as file_handle:
        #up_records = SwissProt.parse(file_handle)

        #        xrefs = defaultdict(list)
        #        xref_ids = list()
        for record in SwissProt.parse(file_handle): #up_records:
            xrefs = defaultdict(list)
            # print(record.entry_name)
            # print(record.data_class)
            # print(record.molecule_type)
            # print(record.accessions)
            # print(record.created)
            # print(record.sequence_update)
            # print(record.annotation_update)
            # print(record.description)
            # print(record.gene_name)
            # print(record.organism)
            # print(record.organelle)
            # print(record.taxonomy_id)
            # print(record.host_taxonomy_id)
            # print(record.references)
            # print(record.comments)
            # print(record.cross_references)
            # print(record.keywords)
            # print(record.features)

            #ref_record = SwissProt._read_rx(record.references,'RX')
            rdf_buffer = ''
            for taxID in record.taxonomy_id:
                if taxID :#in taxon_ids or record.entry_name in lookup_list:
                    # Accession
                    if len(record.accessions) > 1:
                        prim_accession = record.accessions.pop(0)
                        prot_counter += 1
                        rdf_buffer += up_ns + prim_accession + "\n"  # output_writer.write(up_ns + prim_accession + "\n")
                        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n"  # output_writer.write("\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n")
                        rdf_buffer += "\t" + dc_ns + "identifier" + "\t" + '"%s"' % (prim_accession) + " ;\n"
                        rdf_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (prim_accession) + " ;\n"
                        for altID in record.accessions:
                            rdf_buffer += "\t" + owl_ns + "sameAs" + "\t" + up_ns + altID + " ;\n"  # output_writer.write("\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + protein_term + " ;\n")
                            rdf_buffer += "\t" + skos_ns + "altLabel" + "\t" +  '"%s"' % (altID) + " ;\n"  # output_writer.write("\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + protein_term + " ;\n")

                    else:
                        prim_accession = record.accessions[0]
                        prot_counter += 1
                        rdf_buffer += up_ns + prim_accession + "\n"
                        rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n"
                        rdf_buffer += "\t" + dc_ns + "identifier" + "\t" + '"%s"' % (prim_accession) + " ;\n"
                        rdf_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (prim_accession) + " ;\n"

                        # rdf_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + protein_term + " ;\n"

                    # Label
                    rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (record.entry_name) + " ;\n"
                    rdf_buffer += "\t" + skos_ns + "prefLabel" + "\t" + '"%s"' % (record.entry_name) + " ;\n"
                    print("%s\n" % (str(prim_accession)))
                    # Description
                    if record.description:
                        descriptions = record.description.split(';')
                        description = descriptions[0][14:]  # .lstrip('RecName: Full=')
                        rdf_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % (description) + " ;\n"
                    #  Gene Name
                    if record.gene_name:
                        for entry in record.gene_name.split(';'):
                            new_entry = re.sub('\{.+?\}', '', entry)
                            if re.findall("Name=",new_entry):
                                value = new_entry.split('=')[1]
                                for symbol in value.split(','):
                                    symbol = re.sub('\s+', '', symbol)
                                    rdf_buffer += "\t" + skos_ns + "prefSymbol" + "\t" + '"%s"' % (
                                                    symbol) + " ;\n"
                            if re.findall("Synonyms=", new_entry):
                                value = new_entry.split('=')[1]
                                for symbol in value.split(','):
                                    symbol = re.sub('\s+', '', symbol)
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                            if re.findall("OrderedLocusNames=", new_entry):
                                value = new_entry.split('=')[1]
                                for symbol in value.split(','):
                                    symbol = re.sub('\s+', '', symbol)
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                            if re.findall("ORFNames=", new_entry):
                                value = new_entry.split('=')[1]
                                for symbol in value.split(','):
                                    symbol = re.sub('\s+', '', symbol)
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"

                    # Taxon
                    rdf_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxID + " ;\n"
                    # Comments
                    if record.comments:
                        comment_buffer = splitComments(record.comments)
                        raw_comment = ''.join(record.comments)
                        comment = raw_comment.replace('"', '')
                        rdf_buffer += comment_buffer
                        print(comment_buffer)

                    # Keywords
                    #                    print record.keywords
                    if record.keywords:
                        for keyword in record.keywords:
                            keywordURI, cleanKeyword = keyword2URI(keyword)
                            keyword_list.add(cleanKeyword)
                            rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + keywordURI + " ;\n"
                    if record.cross_references:
                        for dbs in record.cross_references:
                            dbname = dbs[0]
                            ids = dbs[1]
                            xrefs[dbname].append(ids)
                    for key in xrefs:
                        if key != "GO":
                            db_namespace = key.lower()
                            for dbid in xrefs[key]:
                                # rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + "<" + up_base_uri + db_namespace + "/" + dbid + ">" + " ;\n"
                                rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + "<" + up_base_uri + db_namespace + "/" + dbid + ">" + " ;\n"
                        if key == "GO":
                            for dbid in xrefs[key]:
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns+ re.sub(':','_',dbid) + " ;\n"

                    # reference citation

                    #print(dir(record.references))
                    for obj in record.references.__iter__():
                        for citation in obj.references:
                            if citation[0] == "PubMed":
                                rdf_buffer += "\t" + dcterms_ns + "references" + "\t" + pubmed_ns+ citation[1]  + " ;\n"
                            if citation[0] == "DOI":
                                rdf_buffer += "\t" + dc_ns + "identifier" + "\t" + '"http://dx.doi.org/%s"' %  (citation[1])   + " ;\n"
                    # Corss references using blank node
                    #                    for key in xrefs:
                    #                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + "[" + "\n"
                    #                        rdf_buffer += "\t" + "\t" +  base_vocab_ns + "dbname" + "\t" + '"%s"' % (key) + " ;\n" #"[" +
                    #                        for dbid in xrefs[key]:
                    #                            rdf_buffer += "\t" + "\t" + base_vocab_ns + "id" + "\t" + '"%s"' % (dbid) + " ;\n"
                    #                        rdf_buffer = re.sub(' ;$', '', rdf_buffer)
                    #                        rdf_buffer += "\t" + "\t" + "]" + " ;\n"

                    rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
                    output_writer.write(rdf_buffer)
                    print("\t".join(keyword_list))
        file_handle.close()
    output_writer.close()
    keyword_writer.write(str(getRDFHeaders()))
    for keyword in keyword_list:
        keywordBuffer = ''
        keywordBuffer = keyword2Triples(keyword,keywordBuffer)
        keyword_writer.write(keywordBuffer)

    print("Number of Proteins: %s\n" % (str(prot_counter)))
    print("*************** UniProt RDF conversion completed ************\n")


#                   pp.pprint(record.cross_references) #taxonomy_id cross_references comments description keywords gene_name molecule_type

# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
#up_dir = '/Users/plarmande/workspace2015/AgroLD.old/test_files/uniprot/uniprot_sample.dat'
#ROOT_DIR = '/Users/plarmande/Downloads'
# upToRDF(up_dir, ROOT_DIR,'/Volumes/LaCie/AGROLD/data_update_2019/Rice_Genome_Hub/uniprot_id.txt')
#upToRDF(up_dir, ROOT_DIR, '/Users/plarmande/workspace2015/datasets/uniprot_id.txt')

# code running on bioinfo-inter
up_dir = sys.argv.pop() # path to the uniprot dataset
ROOT_DIR = sys.argv.pop() # path to the root folder
uniprotid_list = sys.argv.pop() # path
print("%s .... %s ... %s ..." % (up_dir,ROOT_DIR,uniprotid_list))
upToRDF(up_dir,ROOT_DIR,uniprotid_list)
# os.system('cp /scratch/larmande/uniprot.plants.ttl /data3/projects/agrold/uniprot/uniprot.trembl.plants.ttl')