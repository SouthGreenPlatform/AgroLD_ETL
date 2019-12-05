'''
Updated on Dec 4, 2019

@author: larmande
'''
# import sys
# sys.path.append("/home/venkatesan/Downloads/biopython-1.65/Bio")
# TODO chercher dans les synonyms les pattern RAPDB et MSU pour creer une nouvelle relation
# TODO ajouter Prot:uri encodedBy Gene:uri
from Bio import SwissProt

from riceKB.globalVars import base_vocab_ns
from riceKB.globalVars import *
import pprint
import re
import os
from _collections import defaultdict


def upToRDF(up_files, rdf_out_dir, additional_file):  # , output_file

    rdf_file = "uniprot.plants.ttl"
    output_file = os.path.join(rdf_out_dir, rdf_file)
    output_writer = open(output_file, "w")
    rdf_buffer = ''
    prot_counter = 0
    lookup_list = set()
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

    output_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    output_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    output_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    output_writer.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    output_writer.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    output_writer.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    output_writer.write(pr + "\t" + dcterms_ns + "<" + dcterms_uri + "> .\n")
    output_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    output_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    output_writer.write(pr + "\t" + sio_ns + "<" + sio_uri + "> .\n")
    output_writer.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n")
    output_writer.write(pr + "\t" + pubmed_ns + "<" + pubmed_uri + "> .\n")
    output_writer.write(pr + "\t" + up_ns + "<" + uniprot + "> .\n\n")

#    for upfile in up_files:
#        file_handle = open(upfile, "r")
    with open(up_files, 'r') as file_handle:
        up_records = SwissProt.parse(file_handle)

        #        xrefs = defaultdict(list)
        #        xref_ids = list()
        for record in up_records:
            xrefs = defaultdict(list)
            #ref_record = SwissProt._read_rx(record.references,'RX')
            rdf_buffer = ''
            for taxID in record.taxonomy_id:
                if taxID in taxon_ids or record.entry_name in lookup_list:
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

                    # Description
                    if record.description:
                        descriptions = record.description.split(';')
                        description = descriptions[0][14:]  # .lstrip('RecName: Full=')
                        rdf_buffer += "\t" + dcterms_ns + "description" + "\t" + '"%s"' % (description) + " ;\n"
                    #                    print description

                    #  Gene Name
                    #                    print record.gene_name

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
                    #                   taxID

                    # Comments
                    if record.comments:
                        raw_comment = ''.join(record.comments)
                        comment = raw_comment.replace('"', '')
                        rdf_buffer += "\t" + rdfs_ns + "comment" + "\t" + '"%s"' % (comment) + " ;\n"
                    #                   print (comment)

                    # Keywords
                    #                    print record.keywords
                    if record.keywords:
                        for keyword in record.keywords:
                            rdf_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + '"%s"' % (keyword) + " ;\n"
                    # Cross References
                    #                    pp.pprint(record.cross_references[0])

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
                                rdf_buffer += "\t" + dc_ns + "identifiers" + "\t" + pubmed_ns + citation[1] + " ;\n"
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
        file_handle.close()
    output_writer.close()
    print("Number of Proteins: %s\n" % (str(prot_counter)))
    print("*************** UniProt RDF conversion completed ************\n")


#                   pp.pprint(record.cross_references) #taxonomy_id cross_references comments description keywords gene_name molecule_type

# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
# up_dir = '/Volumes/LaCie/AGROLD/argoLD_project_all_data_avril_2017/data/uniport/*.dat'
up_dir = '/Users/plarmande/workspace2015/AgroLD.old/test_files/uniprot/uniprot_sample.dat'
ROOT_DIR = '/Users/plarmande/Downloads'
# upToRDF(up_dir, ROOT_DIR,'/Volumes/LaCie/AGROLD/data_update_2019/Rice_Genome_Hub/uniprot_id.txt')
upToRDF(up_dir, ROOT_DIR, '/Users/plarmande/workspace2015/datasets/uniprot_id.txt')

# code running on bioinfo-inter
# up_dir = ['/scratch/larmande/uniprot_trembl.dat']
# ROOT_DIR='/scratch/larmande'
# upToRDF(up_dir, ROOT_DIR )
# os.system('cp /scratch/larmande/uniprot.plants.ttl /data3/projects/agrold/uniprot/uniprot.trembl.plants.ttl')