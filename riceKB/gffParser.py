from __future__ import with_statement
from globalVars import *
from collections import namedtuple
import gzip
import urllib3
import re
import pprint


__author__ = 'larmande'


#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)
attrib = {}
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        key = re.sub('"','', key)
        value = re.sub('"','', value)
        ret[key] = value
        attrib[key]=value
    return ret

def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    map_ds = list()
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else re.sub('"','', parts[0]),
                "source": None if parts[1] == "." else re.sub('"','', parts[1]),
                "type": None if parts[2] == "." else re.sub('"','', parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else re.sub('"','', parts[6]),
                "phase": None if parts[7] == "." else re.sub('"','', parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            map_ds.append(normalizedInfo)
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            #yield GFFRecord(**normalizedInfo)
        return map_ds


pp = pprint.PrettyPrinter(indent=4)
