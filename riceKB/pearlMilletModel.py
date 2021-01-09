import sys
print (sys.path)
from riceKB.globalVars import *
from riceKB.globalVars import base_vocab_ns
from riceKB.gffParser import *
from riceKB.utils import *
import pprint
import re
import os
import json
import pandas as pd
import numpy as np
import rdflib
from rdflib.graph import Graph

# TODO modifier dans la base les predicats develops_from et has_trait avec developsFrom hasTrait
'''
Created on Jan, 2021
This module is created as part of the AgroLD Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Pearl Millet genome annotation \
datasets to normalize with AgroLD

TODO:
    1) Add documentation
    2) better Error handling

Questions: 
- which URI to define for resources such as genes, transcript and proteins ? \
for now its all after http://www.southgreen.fr/agrold/resource/ 

@author: larmande
'''


__author__  = "larmande"

ROOT_DIR = ''
file_input = ''
prefixes = {}
def milletParser(files, type):
    """
       parse the ensembl rdf file
       :param files: path of the file
       write  objects corresponding to prefix and triples
    """