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

path = '/Users/pierre/Downloads/RO.nt'
#path = '/media/elhassouni/donnees/Noeud-plante-projet/workspace/AgroLD/AgroLD_ETL/test_files/urgi/pseudomolecul_wheat.gff'
g = Graph()
g.parse(path, format="nt")

# qres = g.query(
#     """SELECT *
#        WHERE {
#           ?a ?c ?b .
#        }""")
g.serialize(destination='/Users/pierre/Downloads/RO-new.ttl', format='turtle')

# for row in qres:
#     print("%s %s %s" % row)
#pp.pprint(parseGFF3(path))
