import sys
#print sys.path
from globalVars import *
from globalVars import base_vocab_ns
import pprint
import re
import os
import pandas as pd
import numpy as np
import urllib
import json

'''
Created on Sept, 2019
The ensmbl.rdf module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic functions for handling Esembl RDF datasets to normalize with AgroLD

TODO:
    1) Add documentation
    2) better Error handling
@author: larmande
'''
__author__  = "larmande"




ROOT_DIR = '/Users/plarmande/workspace2015/data/'

ensemb_files = [ROOT_DIR + 'Oryza_sativa_japonica.txt']
ensembl_out = '/Users/plarmande/workspace2015/data/'


pp = pprint.PrettyPrinter(indent=4)

print "***************** Esembl RDF data ********************\n"


print "********************************************************\n\n"

# if __name__ == '__main__':
# 	link = "http://data.gramene.org/v60/genes?q=Os08g0494375"
# 	result = connectionError(link)
# 	print(result.content)