from SPARQLWrapper import SPARQLWrapper, JSON, XML
import sys
from riceKB.globalVars import *
from riceKB.utils import *
from riceKB.globalVars import base_vocab_ns
import pprint
import re
import os

'''
Created on December, 2020
The getAllProperties module is created as part of the AgroLD project.


@author: larmande
'''
ROOT_DIR = sys.argv.pop() # "/Users/pierre/Downloads"
endpoint = "http://agrold.ird.fr:8890/sparql" #"http://agrold.southgreen.fr/sparql" #"http://agrold.ird.fr:8890/sparql"

sparql = SPARQLWrapper(endpoint)
# IRD : http://agrold.ird.fr:8890/sparql
# CIRAD : http://agrold.southgreen.fr/sparql
sparql.setQuery("""
BASE <http://www.southgreen.fr/agrold/>

SELECT distinct ?graph
WHERE {
 GRAPH ?graph {
   ?subject ?predicate ?object.
 }
 filter(REGEX(?graph, "^http://www.southgreen.fr/agrold/"))
}
""")

rdf_file = "getAllProperties.txt"
output_file = os.path.join(ROOT_DIR, rdf_file)
output_writer = open(output_file, "w")
#output_writer.write(str(getRDFHeaders()))

sparql.setReturnFormat(JSON)
results = sparql.query().convert()

for result in results["results"]["bindings"]:
    if "graph" in result:

        graph_name = result["graph"]["value"]
        print(graph_name)
        # sparql.setQuery("""
        #     BASE <http://www.southgreen.fr/agrold/>
        #     PREFIX graph:<""" +
        #                 graph_name +
        #     """>
        #
        #     SELECT distinct ?relation
        #     WHERE {
        #         GRAPH graph: {
        #             ?subject ?relation ?object .
        #         }
        #     }
        #     ORDER BY ?relation
        # """)
        # sparql.setReturnFormat(JSON)
        # results = sparql.query().convert()
        # if "relation" in results["results"]["bindings"]:
        #     ttl_buffer = ''
        #     print("\t".join((result["relation"]["value"], graph_name)) )
        #     ttl_buffer += "\t".join((result["relation"]["value"], graph_name))+"\n"
        #     output_writer.write(ttl_buffer)
output_writer.close()