import sys
from riceKB.globalVars import *
import pprint
import requests
from requests.auth import HTTPDigestAuth
import json
import re
import os
import pandas as pd
import numpy as np
from rdflib import *
from rdflib import Graph, URIRef, Literal
from SPARQLWrapper import *
from pandas.io.json import json_normalize

# from SPARQLWrapper import SPARQLWrapper, JSON, XML, RDFXML

'''
Created on May, 2017
The msuModel module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling MSU data


@author: larmande
'''
__author__  = "larmande"

#TODO 1) Add documentation
#TODO 2) needs code Refactoring
#TODO 3) better Error handling

# Replace with the correct URL
url = "http://data.gramene.org/v60/genes?"
genes = 'q=os11g0559200'

ROOT_DIR='/Users/plarmande/workspace2015/data/'
rdf_writer = open(ROOT_DIR+"result.tsv", "w")
sparql = SPARQLWrapper("http://sparql.southgreen.fr/")
graph = ['http://www.southgreen.fr/agrold/gramene.genes']
sparql.setQuery("""
    PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX graph:<http://www.southgreen.fr/agrold/gramene.genes>
    SELECT distinct ?subject
    WHERE { 
        GRAPH graph: { 
         ?subject ?relation ?object .
        } 
        } 
    """)
# sparql.setReturnFormat(JSON)
# results = sparql.query().convert()
# for result in results["results"]["bindings"]:
#         print('%s' % (result["subject"]["value"]))
#         rdf_writer.write('%s\n' % (result["subject"]["value"]))

# req=requests.get(url+genes)
# print(req.json())
#data= json.dumps(req.json(), indent = 4, sort_keys=True)
#json_normalize(data)
#dataframe= pd.read_json(data)
#print(dataframe)
# data = json.load(req.json())

# df = pd.io.json.json_normalize(req.json())

g = Graph()
ncbi= g.parse('https://rdf.ncbi.nlm.nih.gov/pubmed/22218673')
#result = g.parse("/Users/plarmande/workspace2015/data/temp_graph.ttl",format='turtle')
g.serialize(destination="file:/Users/plarmande/workspace2015/data/temp_graph_validated.ttl", format='turtle')