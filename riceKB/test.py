import riceKB.gffParser
import pprint
from riceKB.globalVars import *
from riceKB.utils import *
import re
import os, sys
import datetime
import json
import urllib
import rdflib
from rdflib.graph import Graph
from rdflib import Graph, plugin
from rdflib.serializer import Serializer

pp = pprint.PrettyPrinter(indent=4)
#handle = open("/Users/pierre/Downloads/RO.ttl", "w")

path = '/Users/pierre/workspace2015/datasets/PlantTFDB'
for filename in os.listdir(path):
    if filename.endswith("_TF_list.ttl"):
        new_filename = "Rdflib_" + filename
        output = os.path.join(path,new_filename)
        path_file= os.path.join(path,filename)
        g = Graph()
        g.parse(path_file, format="turtle")
        g.serialize(destination=output, format='turtle')

#path = '/media/elhassouni/donnees/Noeud-plante-projet/workspace/AgroLD/AgroLD_ETL/test_files/urgi/pseudomolecul_wheat.gff'

# qres = g.query(
#     """SELECT *
#        WHERE {
#           ?a ?c ?b .
#        }""")

# for row in qres:
#     print("%s %s %s" % row)
#pp.pprint(parseGFF3(path))
