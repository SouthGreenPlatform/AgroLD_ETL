'''
Updated on Dec 4, 2019

@author: larmande
'''
# TODO chercher dans les synonyms les pattern RAPDB et MSU pour creer une nouvelle relation
# TODO ajouter Prot:uri encodedBy Gene:uri
from Bio import SwissProt

from globalVars import *
#from globalVars import *
from utils import *
import pprint
import re
import os
import sys
from collections import defaultdict

SPECIAL_CHARACTERS = '[@_!#$%^&*<>?/\\|}{~:.,;-]'
# The pattern needs to be a raw string to handle backslashes properly, and we use re.escape to escape special characters
pattern = '[' + re.escape(SPECIAL_CHARACTERS) + ']'

def keyword2URI(keyword):
    cleanKey=''
    if  re.split('{',keyword):
        cleanKey = re.split('{',keyword)[0]
        cleanKey = re.sub(r'^\s+|\s+$', '', cleanKey)
        cleanKey = re.sub(r'\s', '_', cleanKey)
        cleanKey = re.sub(r'/', '_', cleanKey)
        cleanKey = re.sub(pattern, '_', cleanKey)
    else:
        cleanKey = re.sub(r'^\s+|\s+$','',keyword)
        cleanKey = re.sub(r'\s','_',cleanKey)
        cleanKey = re.sub(pattern, '_', cleanKey)
    # remove _ at the end of the keyword
    if cleanKey.endswith('_'):
        cleanKey = cleanKey[:-1]
    keywordURI =  "<" + base_resource_uri + 'keyword/' + cleanKey + ">"
    return keywordURI, cleanKey

def keyword2Triples(cleanedKey, keywordBuffer):
    # cleanKey = re.sub('^\s+|\s+$', '', keyword)
    keyword = re.sub('_', ' ', cleanedKey)
    keywordBuffer += "<" + base_resource_uri + "keyword/" + cleanedKey + "> \n"
    keywordBuffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Keyword" + " ;\n"
    keywordBuffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (keyword) + " .\n"
    return keywordBuffer

def pubmed2RDF(pubmedid, pubmed_dict,taxID):
    buffer = ''
    buffer += pubmed_ns +  pubmedid + " \n"
    buffer += "\t" + rdf_ns + "type" + "\t" + bibo_ns + "Article" + " ;\n"
    title = pubmed_dict['title']
    title = re.sub(r'\"|\'', '', title)
    buffer += "\t" + dc_ns + "title" + "\t" +  "\""+ title + "\" ;\n"
    for key, value in pubmed_dict['references']:
        if key == 'DOI':
            buffer += "\t" + dcterms_ns + "identifier" + "\t" + "\""+ doi_uri + value + "\" ;\n"
    # buffer += "\t" + skos_ns + "exactMatch" + "\t" + pubmed_ns + "2719677" + " ;\n"
    buffer += "\t" + dc_ns + "creator" + "\t" + "\""+  re.split(',',pubmed_dict['authors'])[0] +  " et al.\"  ;\n"
    year = str(re.findall(r'\.*\(\d+\)\.', pubmed_dict['location'])[0])
    year = re.sub(r'\(|\)|\.', '', year)
    buffer += "\t" + dc_ns + "date" + "\t" + "\""+ year + "\" ;\n"
    buffer += "\t" + prism_ns + "publicationDate" + "\t" + "\"" + year + "\" ;\n"
    buffer += "\t" + dc_ns + "source" + "\t" + "\""+ pubmed_dict['location'] + "\" ;\n"
    publicationName = str(re.split(r'\s\d+', pubmed_dict['location'])[0])
    # volPage = str(re.split('\w+\.*\s', pubmed_dict['location'])[1])
    buffer += "\t" + prism_ns + "publicationName" + "\t" + "\"" + publicationName + "\" ;\n"
    # volume = str(re.split(':', volPage)[0])
    # buffer += "\t" + prism_ns + "volume" + "\t" + volume + " ;\n"
    # buffer += "\t" + up_core_ns + "pages" + "\t" + "name" + " .\n\n"
    buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxID + " ;\n"
    buffer += "\t" + rdf_ns + "seeAlso" +  "\t" + "<" + pubmed_rdf_uri + pubmedid +">. \n\n"
    return buffer
#  a bibo:Article;
#   prism:publicationName "Nature science cell";
#   prism:volume "10";
#   prism:number "11";
#   prism:startingPage "123";
#   prism:endingPage "456";
#   dcterms:date "2015-12-08" ;
#   seeAlso <http://rdf.ncbi.nlm.nih.gov/pubmed/1234567>

def cleanUp(text, title, provenance=False):
    clean_text = text.replace('"', '')
    clean_text = clean_text.replace(title, '')
    clean_text = re.sub(r'^\s+', '', clean_text)
    if provenance:
        clean_text = re.sub(r'{.+?}', '', clean_text)
    return clean_text

def splitComments(comments,accession):
    buffer = ''
    for annotation in comments:
        if 'FUNCTION:' in annotation:
            # (RO_0000085  has function)
            function = cleanUp(annotation,'FUNCTION:')
            buffer += "\t" + obo_ns + "RO_0000085" + "\t" + '"%s"' % (function) + " ;\n"
        elif 'SUBCELLULAR LOCATION:' in annotation:
            # RO_0001025  located in
            location = cleanUp(annotation,'SUBCELLULAR LOCATION:', True)
            if 'Note=' in location:
                location = location.split('Note=')[0]
            location_entries = re.split(r'\.|,|;',location)
            for location_entry in location_entries:
                location_entry = re.sub(r'^\s+|\s+$', '', location_entry)
                if location_entry:
                    buffer += "\t" + obo_ns + "RO_0001025" + "\t" + '"%s"' % (location_entry) + " ;\n"
        elif 'TISSUE SPECIFICITY:' in annotation:
            #  up:Tissue_Specificity_Annotation RO_0002206  expressed in
            tissue = cleanUp(annotation,'TISSUE SPECIFICITY:')
            buffer += "\t" + obo_ns + "RO_0002206" + "\t" + '"%s"' % (tissue) + " ;\n"
        elif 'DOMAIN:'  in annotation:
            #  RO_0002524:contains
            domain = cleanUp(annotation,'DOMAIN:')
            buffer += "\t" + obo_ns + "RO_0002524" + "\t" + '"%s"' % (domain) + " ;\n"
        elif 'SIMILARITY:' in annotation:
            family = cleanUp(annotation,'SIMILARITY:')
            # RO_0002350
            buffer += "\t" + obo_ns + "RO_0002350" + "\t" + '"%s"' % (family) + " ;\n"
        # elif 'SEQUENCE CAUTION:' in annotation:
        #     # Sequence=BAD35485.1; Type=Erroneous gene model prediction; Evidence={ECO:0000305};
        #     # Sequence=BAS98577.1; Type=Erroneous gene model prediction; Evidence={ECO:0000305};
        #     print('found SEQUENCE CAUTION:')
        elif 'CATALYTIC ACTIVITY:' in annotation:
            # hasActivity
            activity = cleanUp(annotation, 'CATALYTIC ACTIVITY:')
            buffer += "\t" + base_vocab_ns + "hasActivity" + "\t" + '"%s"' % (activity) + " ;\n"
        elif 'COFACTOR:'  in annotation:
            cofactors = cleanUp(annotation, 'COFACTOR:')
            # hasCofactor
            cofactor_list = cofactors.split(';')
            for name in cofactor_list:
                if 'Name=' in name:
                    cofactor = re.sub('Name=','',name)
                    cofactor = re.sub(r'\s+','',cofactor)
                    buffer += "\t" + base_vocab_ns + "hasCofactor" + "\t" + '"%s"' % (cofactor) + " ;\n"
        elif 'ACTIVITY REGULATION:'  in annotation:
            regulation = cleanUp(annotation, 'ACTIVITY REGULATION:')
            buffer += "\t" + base_vocab_ns + "activityRegulation" + "\t" + '"%s"' % (regulation) + " ;\n"
            # activityRegulation
        elif 'BIOPHYSICOCHEMICAL PROPERTIES:'  in annotation:
            # bioPhysioChemicalProperties
            biophysic = cleanUp(annotation, 'BIOPHYSICOCHEMICAL PROPERTIES:')
            buffer += "\t" + base_vocab_ns + "biophysiochemicalProperties" + "\t" + '"%s"' % (biophysic) + " ;\n"
        elif 'SUBUNIT:' in annotation:
            #  BFO_0000050
            subunit = cleanUp(annotation, 'SUBUNIT:')
            buffer += "\t" + obo_ns + "BFO_0000051" + "\t" + '"%s"' % (subunit) + " ;\n"
        elif 'PATHWAY:'  in annotation:
            # RO_0000056
            pathway = cleanUp(annotation, 'PATHWAY:')
            buffer += "\t" + obo_ns + "RO_0000056" + "\t" + '"%s"' % (pathway) + " ;\n"
        elif 'DEVELOPMENTAL STAGE:' in annotation:
            # expressedAt
            stage = cleanUp(annotation, 'DEVELOPMENTAL STAGE:')
            buffer += "\t" + base_vocab_ns + "expressedAt" + "\t" + '"%s"' % (stage) + " ;\n"
        elif 'INDUCTION:' in annotation:
            #  RO_0002334  regulated by
            #  RO_0002336 positively regulated by
            induction = cleanUp(annotation, 'INDUCTION:')
            buffer += "\t" + obo_ns + "RO_0002336" + "\t" + '"%s"' % (induction) + " ;\n"
        elif 'PTM:' in annotation:
            pmt = cleanUp(annotation, 'PTM:')
            buffer += "\t" + base_vocab_ns + "isModified" + "\t" + '"%s"' % (pmt) + " ;\n"
        # elif 'RNA EDITING:' in annotation:
        #     # long text
        #     print('found RNA EDITING:')
        # elif 'MASS SPECTROMETRY:' in annotation:
        #     # long text
        #     print('found MASS SPECTROMETRY:')
        elif 'POLYMORPHISM:' in annotation:
            polymorphism = cleanUp(annotation, 'POLYMORPHISM:')
            #  GENO_0000413
            buffer += "\t" + obo_ns + "GENO_0000413" + "\t" + '"%s"' % (polymorphism) + " ;\n"
        elif 'DISEASE:' in annotation:
            disease =  cleanUp(annotation, 'DISEASE:')
            # SIO:000001
            buffer += "\t" + sio_ns + "SIO:000001" + "\t" + '"%s"' % (disease) + " ;\n"
        elif 'DISRUPTION PHENOTYPE:' in annotation:
            disruption = cleanUp(annotation, 'DISRUPTION PHENOTYPE:')
            # SIO:000001
            buffer += "\t" + sio_ns + "SIO:000001" + "\t" + '"%s"' % (disruption) + " ;\n"
        # elif 'ALLERGEN:' in annotation:
        #     # Causes an allergic reaction in human.
        #     print('found ALLERGEN:')
        elif 'INTERACTION:' in annotation:
            interactions = cleanUp(annotation, 'INTERACTION:')
            interaction_list = re.split(r'EBI-\d+;',interactions)
            for intact in interaction_list:
                interactants = re.split(r'NbExp=\d+;',intact)[0]
                second_interactant = interactants.split(':')[0]
                if second_interactant:
                    cleaned_interactant = re.split(';',second_interactant)[1]
                    cleaned_interactant = re.sub(r'\s+','', cleaned_interactant)
                    if cleaned_interactant:
                        if prot_pattern.match(cleaned_interactant):
                            buffer += "\t" + obo_ns + "RO_0002434" + "\t" +  up_ns + cleaned_interactant + " ;\n"
        #         interacts with	RO_0002434
        # elif 'TOXIC DOSE:' in annotation:
        #     # PD(50) is 0.36 mg/kg by injection in blowfly larvae, in
        #     # CC       PubMed:12885226.
        #     print('found TOXIC DOSE:')
        # elif 'BIOTECHNOLOGY:' in annotation:
        #     # text
        #     print('found BIOTECHNOLOGY:')
        # elif 'PHARMACEUTICAL:' in annotation:
        #     # text
        #     print('found PHARMACEUTICAL:')
        # elif 'MISCELLANEOUS:' in annotation:
        #     # text
        #     print('found MISCELLANEOUS:')
        # elif 'SEQUENCE CAUTION:' in annotation:
        #     # text
        #     print('found SEQUENCE CAUTION:')
        # print(annotation)
    return buffer
def upToRDF(up_files, rdf_out_dir,taxon_id,bank_name):  # , output_file

    rdf_file = "uniprot_" +  taxon_id +  "_" + bank_name + ".plants.ttl"
    rdf_keyword_file = "uniprot_" +  taxon_id +  "_" + bank_name  + ".keyword.ttl"
    rdf_pubmed_file = "uniprot_" +  taxon_id +  "_" + bank_name  + ".pubmed.ttl"
    prot_gene_mapping = "mappings_" +  taxon_id + "_" + bank_name +  "_gene_prot.tsv"
    output_file = os.path.join(rdf_out_dir, rdf_file)
    output_key_file = os.path.join(rdf_out_dir, rdf_keyword_file)
    output_pub_file = os.path.join(rdf_out_dir, rdf_pubmed_file)
    prot_gene_mapping_file = os.path.join(rdf_out_dir, prot_gene_mapping)
    output_writer = open(output_file, "w")
    keyword_writer = open(output_key_file, "w")
    pubmed_writer = open(output_pub_file, "w")
    prot_gene_writer = open(prot_gene_mapping_file, 'w')
    keywordBuffer = ''
    rdf_buffer = ''
    gene_buffer = ''
    pubmed_dict = {}
    prot_counter = 0
    lookup_list = set()
    keyword_list = set()
    pp = pprint.PrettyPrinter(indent=4)
    up_base_uri = "http://purl.uniprot.org/"
    #    up_base_ns = "uniprot_base:"
    # if (additional_file):
        # if (os.path.isfile(str(additional_file))):
        #     with open(additional_file, 'r') as infile:
        #         for prot in infile:
        #             prot = re.sub(r'\s+', '', prot)
        #             lookup_list.add(prot)

    print("************* Converting Uniprot data to RDF ***************\n")

    output_writer.write(str(getRDFHeaders()))
#    for upfile in up_files:
#        file_handle = open(upfile, "r")
    with open(up_files, 'r') as file_handle:
        #up_records = SwissProt.parse(file_handle)

        #        xrefs = defaultdict(list)
        #        xref_ids = list()
            # for record in SeqIO.parse(file_handle, "uniprot-xml"):
            #     print(record.id)
        for record in SwissProt.parse(file_handle): #up_records:

            xrefs = defaultdict(list)
            #ref_record = SwissProt._read_rx(record.references,'RX')
            rdf_buffer = ''
            prot_gene_buffer = ''
            for taxID in record.taxonomy_id:
                if taxID in taxon_ids:
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
                    #print("%s\n" % (str(prim_accession)))
                    # Description
                    #
                    if record.description:
                        descriptions = record.description.split(';')
                        full_name = descriptions.pop(0)
                        full_name = full_name[14:]  # .lstrip('RecName: Full=')
                        full_name = re.sub(r'{.+?}', '', full_name)
                        full_name = full_name.strip()
                        rdf_buffer += "\t" + skos_ns + "prefLabel" + "\t" + '"%s"' % (full_name) + " ;\n"
                        rdf_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % (full_name) + " ;\n"
                        for altName in descriptions:
                            altName = altName.strip()
                            if altName.startswith('AltName: Full='):
                                altName = altName[14:]
                                altName = re.sub(r'{.+?}', '', altName)
                                altName = altName.strip()
                                rdf_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (altName) + " ;\n"
                            if altName.startswith('Short='):
                                altName = altName[6:]
                                altName = re.sub(r'{.+?}', '', altName)
                                altName = altName.strip()
                                rdf_buffer += "\t" + skos_ns + "altLabel" + "\t" + '"%s"' % (altName) + " ;\n"
                        #.lstrip('AltName: Full=')

                    #  Gene Name
                    if record.gene_name:
                        for entry_dic in record.gene_name:
                            if entry_dic.get('Name'):
                                gene_name = entry_dic['Name']
                                new_entry = re.sub(r'{.+?}', '', gene_name)
                                new_entry = re.sub(r'\s+', '', new_entry)
                                new_entry = re.sub(r'\"+', '', new_entry)
                                new_entry = new_entry.strip()
                                rdf_buffer += "\t" + skos_ns + "prefSymbol" + "\t" + '"%s"' % (
                                                    new_entry) + " ;\n"
                            if entry_dic.get('Synonyms'):
                                for symbol in entry_dic['Synonyms']:
                                    symbol = re.sub(r'{.+?}', '', symbol)
                                    symbol = re.sub(r'\s+', '', symbol)
                                    symbol = re.sub(r'\"+', '', symbol)
                                    symbol = symbol.strip()
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                            if entry_dic.get('OrderedLocusNames'):
                                # remove the curly braces and the content within them
                                OrderedLocusNames = re.sub(r'{.+?}', '', str(entry_dic['OrderedLocusNames']))
                                # convert the string to a list
                                OrderedLocusNames = OrderedLocusNames.split(',')
                                # loop through the list
                                for symbol in OrderedLocusNames:
                                    symbol = re.sub(r'\[', '', symbol)  # remove the square brackets
                                    symbol = re.sub(r'\]', '', symbol)  # remove the square brackets
                                    symbol = re.sub(r'\'', '', symbol)  # remove the cotes
                                    symbol = re.sub(r'{.+?}', '', symbol)
                                    symbol = re.sub(r'\s+', '', symbol)
                                    symbol = re.sub(r'\"+', '', symbol)
                                    symbol = symbol.strip()
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                                    rdf_buffer += "\t" + base_vocab_ns + "orderedLocusName" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                                    gene_buffer +=  prim_accession + "\t"  +  symbol + "\n"
                            if entry_dic.get('ORFNames'):
                                for symbol in entry_dic['ORFNames']:
                                    symbol = re.sub(r'{.+?}', '', symbol)
                                    symbol = re.sub(r'\s+', '', symbol)
                                    symbol = re.sub(r'\"+', '', symbol)
                                    symbol = symbol.strip()
                                    rdf_buffer += "\t" + skos_ns + "altSymbol" + "\t" + '"%s"' % (
                                        symbol) + " ;\n"
                    # Taxon
                    rdf_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxID + " ;\n"
                    # Comments
                    if record.comments:
                        comment_buffer = splitComments(record.comments,prim_accession)
                        raw_comment = ''.join(record.comments)
                        comment = raw_comment.replace('"', '')
                        rdf_buffer += comment_buffer
                        # print(comment_buffer)

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
                        if key != "GO" and key != "InterPro":
                            db_namespace = key.lower()
                            for dbid in xrefs[key]:
                                # rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + "<" + up_base_uri + db_namespace + "/" + dbid + ">" + " ;\n"
                                # clean the dbid from space characters
                                clean_dbid = re.sub(r'\s+', '', dbid)
                                #clean_dbid = re.sub(pattern, '_', dbid)
                                #rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + "<" + up_base_uri + db_namespace + "/" + clean_dbid + ">" + " ;\n"
                        elif key == "GO":
                            for dbid in xrefs[key]:
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns+ re.sub(':','_',dbid) + " ;\n"
                        elif key == "InterPro":
                            for dbid in xrefs[key]:
                                rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + interpro_ns+ re.sub(':','_',dbid) + " ;\n"

                    # reference citation

                    # print(dir(record.references))
                    for obj in record.references.__iter__():
                        for citation in obj.references:
                            if citation[0] == "PubMed":
                                pubmed_dict[citation[1]] = {'title': obj.title, 'authors': obj.authors, 'location': obj.location, 'references':obj.references}
                                rdf_buffer += "\t" + dcterms_ns + "references" + "\t" + pubmed_ns+ citation[1]  + " ;\n"
                            # if citation[0] == "DOI":
                            #     rdf_buffer += "\t" + dc_ns + "identifier" + "\t" + '"http://dx.doi.org/%s"' %  (citation[1])   + " ;\n"
                    # Corss references using blank node
                    #                    for key in xrefs:
                    #                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + "[" + "\n"
                    #                        rdf_buffer += "\t" + "\t" +  base_vocab_ns + "dbname" + "\t" + '"%s"' % (key) + " ;\n" #"[" +
                    #                        for dbid in xrefs[key]:
                    #                            rdf_buffer += "\t" + "\t" + base_vocab_ns + "id" + "\t" + '"%s"' % (dbid) + " ;\n"
                    #                        rdf_buffer = re.sub(' ;$', '', rdf_buffer)
                    #                        rdf_buffer += "\t" + "\t" + "]" + " ;\n"

                    rdf_buffer = re.sub(' ;$', ' .\n', rdf_buffer)
                    RDF_validation(rdf_buffer,output_writer,prim_accession,ROOT_DIR)
                    #output_writer.write(rdf_buffer)
        file_handle.close()
    output_writer.close()
    keyword_writer.write(str(getRDFHeaders()))
    pubmed_writer.write(str(getRDFHeaders()))
    prot_gene_writer.write(gene_buffer)
    for keyword in keyword_list:
        keywordBuffer = ''
        keywordBuffer = keyword2Triples(keyword,keywordBuffer)
        RDF_validation(keywordBuffer, keyword_writer, keyword, ROOT_DIR)
        #keyword_writer.write(keywordBuffer)
    for pubmedid in pubmed_dict.keys():
        pubmedBuffer = ''
        pubmedBuffer = pubmed2RDF(pubmedid, pubmed_dict[pubmedid],taxon_id)
        RDF_validation(pubmedBuffer, pubmed_writer, pubmedid, ROOT_DIR)
        #pubmed_writer.write(pubmedBuffer)
    keyword_writer.close()
    pubmed_writer.close()
    print("Number of Proteins: %s\n" % (str(prot_counter)))
    print("*************** UniProt RDF conversion completed ************\n")


#                   pp.pprint(record.cross_references) #taxonomy_id cross_references comments description keywords gene_name molecule_type

# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
#up_dir = '/Users/plarmande/workspace2015/AgroLD.old/test_files/uniprot/uniprot_sample.dat'
#ROOT_DIR = '/Users/plarmande/Downloads'
# upToRDF(up_dir, ROOT_DIR,'/Volumes/LaCie/AGROLD/data_update_2019/Rice_Genome_Hub/uniprot_id.txt')
#upToRDF(up_dir, ROOT_DIR, '/Users/plarmande/workspace2015/datasets/uniprot_id.txt')
# create a function that take rdf_buffer as input and validate RDF with RDFLib



# code running on bioinfo-inter
up_dir = sys.argv.pop() # path to the uniprot dataset
ROOT_DIR = sys.argv.pop() # path to the root folder
taxon_id = sys.argv.pop() # taxon id
bank_name = sys.argv.pop() # bank name
print("%s .... %s ... %s ....%s" % (up_dir,ROOT_DIR,taxon_id,bank_name))
upToRDF(up_dir,ROOT_DIR,taxon_id,bank_name)
# os.system('cp /scratch/larmande/uniprot.plants.ttl /data3/projects/agrold/uniprot/uniprot.trembl.plants.ttl')