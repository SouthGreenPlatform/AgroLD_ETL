from __future__ import with_statement
from riceKB.globalVars import *
from collections import namedtuple
import gzip
import urllib3
import re
import pprint


__author__ = 'larmande'


#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)
attrib = {}
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        key = re.sub('"','', key)
        value = re.sub('"','', value)
        ret[key] = value
        attrib[key]=value
    return ret

def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    map_ds = list()
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else re.sub('"','', parts[0]),
                "source": None if parts[1] == "." else re.sub('"','', parts[1]),
                "type": None if parts[2] == "." else re.sub('"','', parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else re.sub('"','', parts[6]),
                "phase": None if parts[7] == "." else re.sub('"','', parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            map_ds.append(normalizedInfo)
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            #yield GFFRecord(**normalizedInfo)
        return map_ds


def RDFConverter(ds, output_file):
    os_japonica_buffer = ''  # initilised the buffer at zero
    number_match_part_sbgi = 0
    number_match_part_kome = 0
    number_exon = 0
    number_cds = 0
    line_number = 0
    number_five_prime_UTR = 0
    number_three_prime_UTR = 0
    rdf_writer = open(output_file, "w")
    chromosome_list = list()
    gene_list = list()
    mRNA_list = list()
    taxon_id = "39947"
    source_project = "JGI"
    schema_number = "3.1"
    ssp = "kitaake/"
    print("************* RDF conversion begins***********\n")
    rdf_writer.write(base + "\t" + "<" + base_uri + "> .\n")
    rdf_writer.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    rdf_writer.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    rdf_writer.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    rdf_writer.write(pr + "\t" + owl_ns + "<" + owl_uri + "> .\n")
    rdf_writer.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n")
    rdf_writer.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    rdf_writer.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    rdf_writer.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    rdf_writer.write(pr + "\t" + protein_ns + "<" + protein_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    rdf_writer.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n")
    rdf_writer.write(pr + "\t" + interpro_ns + "<" + interpro_uri + "> .\n")
    rdf_writer.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n")
    rdf_writer.write(pr + "\t" + "Phytozome:" + "<" + identifiers_uri + "phytosome/> .\n")
    rdf_writer.write(pr + "\t" + "PFAM:" + "<" + identifiers_uri + "pfam/> .\n")
    rdf_writer.write(pr + "\t" + "Panther:" + "<" + identifiers_uri + "panther.family/> .\n")
    rdf_writer.write(pr + "\t" + "KOG:" + "<" + identifiers_uri + "kog/> .\n")
    rdf_writer.write(pr + "\t" + "TAIR:" + "<" + 	"http://arabidopsis.org/servlets/TairObject?accession=" + "> .\n")
    rdf_writer.write(pr + "\t" + "KEGG:" + "<" + identifiers_uri + "kegg/> .\n")
    rdf_writer.write(pr + "\t" + "EC:" + "<" + identifiers_uri + "ec-code/> .\n")
    rdf_writer.write(pr + "\t" + "InterPro:" + "<" + identifiers_uri + "interpro/> .\n")
    rdf_writer.write(pr + "\t" + "TrEMBL:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + "MSU:" + "<" + identifiers_uri + "ricegap/> .\n")
    rdf_writer.write(pr + "\t" + "SwissProt:" + "<" + identifiers_uri + "uniprot/> .\n")
    rdf_writer.write(pr + "\t" + faldo_ns + "<" + faldo + "> .\n")
    rdf_writer.write(pr + "\t" + base_resource_ns + "<" + base_resource_uri + "> .\n\n")

    for records in ds:
        line_number += 1
        # loop over the file
        #if records['type'] == "chromosome":
        if not records['seqid'] in chromosome_list:
            genome_buffer = ""
            chromosome_list.append(records['seqid'])
            chromosome_nb = getChromosomeNumber(records['seqid'])
            if 'Chr' in records['seqid']:
                genome_buffer += "<" + chromosome_uri + taxon_id +"/"+ ssp + chromosome_nb+ ">\n"
                genome_buffer += "\t" +  obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inAssembly" + "\t" + "\"Reference-" + source_project + "\" ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inSchemaNumber" + "\t" + "\""+ schema_number + "\"" + " ;\n"
                genome_buffer = re.sub(' ;$', ' .\n\n', genome_buffer)
                rdf_writer.write(genome_buffer)
            else:
                genome_buffer += "<" + chromosome_uri + taxon_id + "/" + ssp + chromosome_nb + ">\n"
                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + obo_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inAssembly" + "\t" + "\"Reference-" + source_project + "\" ;\n"
                genome_buffer += "\t" + base_vocab_ns + "inSchemaNumber" + "\t" + "\"" + schema_number + "\"" + " ;\n"
                genome_buffer = re.sub(' ;$', ' .\n\n', genome_buffer)
                rdf_writer.write(genome_buffer)

        # filtering for gene entries
        if records['type'] == "gene":

            if not records['attributes']['Name'] in gene_list:
                genome_buffer = ''
                (strand,position) = getStrandValue(records['strand'])
                strand = str(strand)
                # print the corresponding gene associated at mRNAs
                # os_japonica_buffer = ''
                gene_list.append(records['attributes']['Name'])
                genome_buffer += base_resource_ns + records['attributes']['Name'] + "\n"
                genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                if 'Note' in records['attributes']:
                    genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes']['Note'] + "\" ;\n"
                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                                 + chromosome_nb + ":" + \
                                 str(records['start']) + "-" + str(records['end']) + ":" + \
                                 strand + "> ;\n"
                genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
                genome_buffer += getFaldoRegion(taxon_id, ssp, chromosome_nb, records['start'], records['end'],
                                                records['strand'])
                rdf_writer.write(genome_buffer)
        # filtering for mRNA
        if records['type'] == "mRNA":

            if not records['attributes']['ID'] in mRNA_list:
                genome_buffer = ''
                (strand,position) = getStrandValue(records['strand'])
                strand = str(strand)
                go_list = list()
                # print the corresponding gene associated at mRNAs
                mRNA_list.append(records['attributes']['ID'])
                genome_buffer += base_resource_ns + records['attributes']['ID'] + "\n"
                genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
                genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "mRNA" + " ;\n"
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
                if 'Note' in records['attributes']:
                    genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                        'Note'] + "\" ;\n"

                genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
                if 'Parent' in records['attributes']:
                    parent = records['attributes']['Parent'].split('.')[0]
                    genome_buffer += "\t" + base_vocab_ns + "developsFrom" + "\t" + base_resource_ns + parent +";\n"
                if 'Dbxref' in records['attributes']:
                    for terms in records['attributes']['Dbxref'].split(','):
                        genome_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" +  terms + " ;\n"
                if 'Ontology_term' in records['attributes']:
                    for terms in records['attributes']['Ontology_term'].split(','):
                        if 'GO:' not in terms:
                            genome_buffer += "\t" + base_vocab_ns + "hasAnnotation" + "\t" + "\"" + terms.split(':')[1] + "\" ;\n"
                        else:
                            if terms not in go_list:
                                genome_buffer += "\t" + base_vocab_ns + "classifiedWith" + "\t" + obo_ns + \
                                             re.sub(':', '_', terms) + " ;\n"
                                go_list.append(terms)
                genome_buffer += "\t" + faldo_ns + "location" + "\t\t" + "<" + chromosome_uri + taxon_id +"/"+ ssp \
                                 + chromosome_nb + ":" + \
                                 str(records['start']) + "-" + str(records['end']) + ":" + \
                                 strand + "> ;\n"
                genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
                genome_buffer += getFaldoRegion(taxon_id, ssp, chromosome_nb,records['start'],records['end'],records['strand'])
                rdf_writer.write(genome_buffer)

        if records['type'] == "polypeptide":
            genome_buffer = ''
            # print the corresponding gene associated at mRNAs
            # os_japonica_buffer = ''
            gene_list.append(records['attributes']['ID'])
            genome_buffer += base_resource_ns + records['attributes']['ID'] + "\n"
            genome_buffer += "\t" + base_vocab_ns + "sourceProject" + "\t" + " \"" + source_project + "\" ;\n"
            genome_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Protein" + " ;\n"
            if 'Name' in records['attributes']:
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['Name'] + "\" ;\n"
            else:
                genome_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + records['attributes']['ID'] + "\" ;\n"
            if 'Note' in records['attributes']:
                genome_buffer += "\t" + dcterms_ns + "description" + "\t" + " \"" + records['attributes'][
                    'Note'] + "\" ;\n"
            if 'Derives_from' in records['attributes']:
                genome_buffer += "\t" + base_vocab_ns + "derivesFrom" + "\t"  + base_resource_ns + records['attributes']['Derives_from'] + " ;\n"
            genome_buffer += "\t" + obo_ns + "RO_0002162" + "\t\t" + ncbi_tax_ns + taxon_id + " ;\n"
            genome_buffer = re.sub(' ;$', ' .\n', genome_buffer)
            rdf_writer.write(genome_buffer)
        # break

print(attrib)
pp = pprint.PrettyPrinter(indent=4)
path = '/Users/pierre/workspace2015/datasets/Oryza_sativa_Kitaake_3.1.gff3'
path_output = '/Users/pierre/workspace2015/datasets/Oryza_sativa_Kitaake_3.1.ttl'
ds= parseGFF3(path)

RDFConverter(ds, path_output)