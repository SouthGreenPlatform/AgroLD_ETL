import riceKB.gffParser
import pprint
import pprint
from riceKB.globalVars import *
from riceKB.utils import *
import re
import os, sys
import datetime
import json
import pandas as pd
import numpy as np
import rdflib
from rdflib.graph import Graph

pp = pprint.PrettyPrinter(indent=4)



path = '/Users/plarmande/workspace2015/datasets/OryzabaseGeneListEn_20200114010108.ttl'
#path = '/media/elhassouni/donnees/Noeud-plante-projet/workspace/AgroLD/AgroLD_ETL/test_files/urgi/pseudomolecul_wheat.gff'
g = Graph()
g.parse(path, format="turtle")

qres = g.query(
    """SELECT *
       WHERE {
          ?a ?c ?b .
       }""")

# for row in qres:
#     print("%s %s %s" % row)
#pp.pprint(parseGFF3(path))
