#!/usr/bin/env python
'''
Created on June 8, 2020
The stringDB module is created as part of the AgroLD project.

This module contains Parsers, RDF converters and generic functions for handling StringDB data

#TODO:
    1) Add documentation
    2) Fix Gramene record trailing space in the parser, now it is being handled in the RDF converter
    3) better Error handling
@author: venkatesan

'''
# header of taxonid.protein.links.v.txt
# protein1 protein2 combined_score
# 4530.OS01T0100100-01 4530.OS05T0466900-01 208
# header of taxonid.protein.actions.v.txt

# item_id_a	item_id_b	mode	action	is_directional	a_is_acting	score
# 4530.OS01T0100100-01	4530.OS01T0179700-01	activation	activation	t	t	227
# 4530.OS01T0100100-01	4530.OS01T0179700-01	binding		f	f	194
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
