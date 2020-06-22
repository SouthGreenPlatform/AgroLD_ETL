from SPARQLWrapper import SPARQLWrapper, JSON, XML
import sys
from riceKB.globalVars import *
from riceKB.utils import *
from riceKB.globalVars import base_vocab_ns
import pprint
import re
import os

'''
Created on June, 2020
The uniprotEnrichment module is created as part of the AgroLD project.


@author: larmande
'''
ROOT_DIR = "/Users/plarmande/Downloads"
sparql = SPARQLWrapper("http://agrold.ird.fr:8890/sparql")
# IRD : http://agrold.ird.fr:8890/sparql
# CIRAD : http://agrold.southgreen.fr/sparql
sparql.setQuery("""
BASE <http://www.southgreen.fr/agrold/>
PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>
PREFIX vocab: <vocabulary/>
PREFIX resource:<resource/>
PREFIX obo:<http://purl.obolibrary.org/obo/>
PREFIX uniprot:<http://purl.uniprot.org/uniprot/>
PREFIX sio:<http://semanticscience.org/resource/>


SELECT ?protein_id ?label ?symbol FROM <http://www.southgreen.fr/agrold/uniprot.plants>
WHERE {              
?protein_id ?p ?o;       
rdfs:label ?label ;
obo:RO_0002162  <http://identifiers.org/taxonomy/39947>  .
optional {       
?protein_id skos:altSymbol ?symbol.
}
}

""")

rdf_file = "uniprotEnriched.ttl"
output_file = os.path.join(ROOT_DIR, rdf_file)
output_writer = open(output_file, "w")
output_writer.write(str(getRDFHeaders()))

sparql.setReturnFormat(JSON)
results = sparql.query().convert()
for result in results["results"]["bindings"]:

    gene_uri = "http://www.southgreen.fr/agrold/resource/" + result["symbol"]["value"]
    sparql.setQuery("""
        ASK WHERE { <""" +
                    gene_uri +
                    """> 
                    rdf:type <http://www.southgreen.fr/agrold/vocabulary/Gene>
        }    
    """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    if results["boolean"]:
        ttl_buffer = ''
        #print("\t".join((result["protein_id"]["value"],result["label"]["value"], result["symbol"]["value"])) )
        ttl_buffer += "<" + result["protein_id"]["value"] + ">\t" + sio_ns + "SIO_000339" + "\t<" + gene_uri + ">.\n"
        output_writer.write(ttl_buffer)

output_writer.close()