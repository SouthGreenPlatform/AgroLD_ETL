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
Created on June, 2021
This module is created as part of the AgroLD Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Ensembl genome annotation \
datasets to normalize with AgroLD

It runs with several files downloaded from Ensembl Plants project available at  ftp.ensemblgenomes.org


@author: larmande
@author: larmande
'''


__author__  = "larmande"





#TEST PARAM
path = '/Users/plarmande/workspace2015/datasets/Oryza_sativa_Kitaake_3.1.gff3'
path_output = '/Users/plarmande/workspace2015/datasets/Oryza_sativa_Kitaake_3.1.ttl' # The output
#path = '/opt/TOS_DI-20141207_1530-V5.6.1/workspace/gff_data_orygeneDB/os_japonica/os_indicaCancat.gff3'    # The input
#path_output = '/home/elhassouni/Bureau/japonica.ttl' # The output
ds = parseGFF3(path)   # The parsing file
pp.pprint(ds)    # For to see in teminal the parsing

#os_indicaModele(ds, path_output)  # The path_output)  # The tranformation fonction tropGeneToRdf(input, output)

RDFConverter(ds, path_output)
