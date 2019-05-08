#!/usr/bin/env python
'''
Created on Jun 19, 2014
Updated on Apr 17, 2019
The grameneParsers module is created as part of the Rice Knowledge Base project.

This module contains Parsers, RDF converters and generic tools for handling Gramene data

TODO : fix taxon assignation currently based on file naming
TODO : use panda DF to parse tsv file and use more generic headers to build the dictionnary
T
    1) Add documentation
    2) better Error handling
@author: larmande
'''
import pprint
from globalVars import *
from globalVars import base_vocab_ns
import re
import os
import pandas as pd
import numpy as np
#from __builtin__ import map

'''
 Parsers 
''' 
'''
def genomeParser(infile):
    pp = pprint.PrettyPrinter(indent=4)
    geneHash = {}

    fileReader = open(infile, "r")
    lines = fileReader.readlines()
    lines.pop(0) # remove header
#    pp.pprint(lines[0])
    for line in lines:
        prot_ids = list()
        line = re.sub('\n$', '', line)
        records = line.split('\t')
#        print records[0]
        #pp.pprint(records[0])
        if len(records) == 10:
            geneId = records.pop(0)
#            prot_id = list()
            if geneId not in geneHash:
                geneHash[geneId] = {'Name':records[0],
                                    'Description':records[1],
                                    'Chromosome': records[2],
                                    'Start': records[3],
                                    'End': records[4],
                                    'Biotype': records[5],
                                    'ProtID': [],
                                    'GO': {}
                                    }
            if geneId in geneHash:
                if records[6]:
                    prot_ids.append(records[6])                                  
                geneHash[geneId]['ProtID'].extend(prot_ids)  #records[6]
                geneHash[geneId]['GO'][records[7]] = records[8]  

        if len(records) == 11:
            geneId = records.pop(0)
            if geneId not in geneHash:
                geneHash[geneId] = {'Name':records[0],
                                    'Description':records[1],
                                    'Chromosome': records[2],
                                    'Start': records[3],
                                    'End': records[4],
                                    'Biotype': records[5],                                    
                                    'GeneID': records[6],
                                    'ProtID': [],
                                    'GO': {}
                                    }

            if geneId in geneHash:
                if records[7]:
                    prot_ids.append(records[7])
                geneHash[geneId]['ProtID'].extend(prot_ids) #records[6]
                geneHash[geneId]['GO'][records[8]] = records[9]  

        if len(records) == 12:
            geneId = records.pop(0)
            loci = list()
            if geneId not in geneHash:
                geneHash[geneId] = {'Name':records[0],
                                    'Description':records[1],
                                    'Chromosome': records[2],
                                    'Start': records[3],
                                    'End': records[4],
                                    'Biotype': records[5],                                    
                                    'GeneID': records[6],
                                    'TAIRlocus': '',
                                    'RapID': '',
                                    'TIGRlocus': [],                                     
                                    'ProtID': [],
                                    'GO': {}
                                    }
                    
            if geneId in geneHash:
                if records[7]:
                    tigr_loci_pattern = re.match(r'^LOC_', records[7]) #(^LOC_)|(^No)
                    rap_pattern = re.match(r'^Os', records[7]) 
                    tair_pattern = re.match(r'^AT', records[7])     
                    if tigr_loci_pattern:
                        loci.append(records[7])
                        geneHash[geneId]['TIGRlocus'].extend(loci)
                    if rap_pattern:
                        geneHash[geneId]['RapID'] = records[7]
                    if tair_pattern:
                        geneHash[geneId]['TAIRlocus'] = records[7]
                if records[8]:
                    prot_ids.append(records[8])
                    geneHash[geneId]['ProtID'].extend(prot_ids) #records[6]
                if records[9] and records[10]:          
                    geneHash[geneId]['GO'][records[9]] = records[10]     # for A.thaliana dataset records[9],records[10] contains PO terms and domain respectively  
    fileReader.close()
    return geneHash
  '''
def geneParser(infile):
    
#    pp = pprint.PrettyPrinter(indent=4)
    gene_hash = {}
    tigr_pattern = re.compile(r'^LOC\_Os\d{1,2}g\d{5}\.\d$')
    rap_pattern = re.compile(r'^Os\d{2}g\d{7}$')
    tair_pattern = re.compile(r'^AT[1-5]G\d{5}$')
    prot_pattern = re.compile(r'^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')
    ont_pattern = re.compile(r'^\w+\:\d{7}$')
    
    file_reader = open(infile, "r")
    lines = file_reader.readlines()
    lines.pop(0) # remove header
    
    for line in lines:
        line = re.sub('\n$', '', line)
        records = line.split('\t')

        gene_id = records.pop(0)
        print(len(records))
        # Building data structure
        if gene_id not in gene_hash:
            gene_hash[gene_id] = {
                                  'Name': records[0],
                                  'Description': records[1],
                                  'Chromosome': records[3],
                                  'Start': records[4],
                                  'End': records[5],
                                  'Strand': records[6],
                                  'Biotype': records[2],
                                  'RapID': '',
                                  'TairLocus': '',
                                  'TIGRlocus': {},
                                  'ProtID': {},
                                  'Ontology': {},
                                  'StringID':''
                                  }

        if gene_id in gene_hash:
            # Records of sps:  O.s.japonica
            if len(records) == 12 and rap_pattern.match(gene_id):
                if rap_pattern.match(gene_id):
                    gene_hash[gene_id]['RapID'] = gene_id
                if records[7]:
                    if prot_pattern.match(records[7]):
                        gene_hash[gene_id]['ProtID'][records[7]] = '-'
                if records[8]:
                    if prot_pattern.match(records[8]):
                        gene_hash[gene_id]['ProtID'][records[8]] = '-'
                if records[9]:
                    gene_hash[gene_id]['Ontology'][records[10]] = records[9]
                if records[11]:
                    gene_hash[gene_id]['StringID'] = records[11]
            else: # Arabidopsis
                if tair_pattern.match(gene_id):
                    if records[7]:
                        gene_hash[gene_id]['TairLocus'] = records[7]
                    if records[8]:
                        gene_hash[gene_id]['ProtID'][records[8]] = '-'
                    if records[9]:
                        gene_hash[gene_id]['ProtID'][records[9]] = '-'
                    if records[10]:
                        gene_hash[gene_id]['Ontology'][records[11]] = records[10]

            # Records of sps: & O. Glaberrima & O.s.indica, T. Aestivum, S.bicolor, T. urartu, O. Nivara
            if len(records) == 11:

                if records[7]:
                    gene_hash[gene_id]['ProtID'][records[7]] = '-'
                if records[8]:
                    gene_hash[gene_id]['ProtID'][records[8]] = '-'
                if records[9]:
                    if ont_pattern.match(records[9]):
                        gene_hash[gene_id]['Ontology'][records[10]] = records[9]
                    else:
                        gene_hash[gene_id]['ProtID'][records[9]] = '-'
                    #                        prot_list.append(records[7])
                    # gene_hash[gene_id]['Ontology'][records[10]] = records[9]
            # Records of sps: O. Glumaepatula, Oryza_punctata
            if len(records) == 10:
                if records[7]:
                    if prot_pattern.match(records[7]):
                        gene_hash[gene_id]['ProtID'][records[7]] = '-'
                if records[8]:
                    if prot_pattern.match(records[8]):
                        gene_hash[gene_id]['ProtID'][records[8]] = '-'
                    if ont_pattern.match(records[8]):
                        gene_hash[gene_id]['Ontology'][records[9]] = records[8]
                if records[9]:
                    if ont_pattern.match(records[9]):
                        gene_hash[gene_id]['Ontology'][records[10]] = records[9]
            # Records of sps: O.barthii, O.meridionalis, O.longistaminata,
            if len(records) == 9:
                if records[6]:
                    if tigr_pattern.match(records[6]):
                        gene_hash[gene_id]['TIGRlocus'][records[6]] = '-'
#                        tigr_loci.append(records[6])
#                        gene_hash[gene_id]['TIGRlocus'].extend(tigr_loci)
                    if rap_pattern.match(records[6]):
                        gene_hash[gene_id]['RapID'] = records[6]
#                    elif tair_pattern.match(records[6]):
#                        gene_hash[gene_id]['Tairlocus'] = records[6]
                    if prot_pattern.match(records[6]):
                        gene_hash[gene_id]['ProtID'][records[6]] = '-'
#                        prot_list.append(records[6])
#                        gene_hash[gene_id]['ProtID'].extend(prot_list)
                if records[7]:
                    if ont_pattern.match(records[7]):
                        gene_hash[gene_id]['Ontology'][records[8]] = records[7]
                    else:
                        gene_hash[gene_id]['ProtID'][records[7]] = '-'
#                        prot_list.append(records[7])
#                        gene_hash[gene_id]['ProtID'].extend(prot_list)
                if records[8]:
                    if prot_pattern.match(records[8]):
                        gene_hash[gene_id]['ProtID'][records[8]] = '-'
#                        prot_list.append(records[8])
#                        gene_hash[gene_id]['ProtID'].extend(prot_list)
            # Records of sps: O.glaberrima, T.aestivum, T.urartu, Z.mays 
            if len(records) == 8:
                if records[5]:
                    gene_hash[gene_id]['ProtID'][records[5]] = '-'
#                    prot_list.append(records[6])
#                    gene_hash[gene_id]['ProtID'].extend(prot_list)
                if records[6]:
                    gene_hash[gene_id]['Ontology'][records[6]] = records[7]
                    #gene_hash[gene_id]['ProtID'][records[7]] = '-'
#                    prot_list.append(records[7])
#                    gene_hash[gene_id]['ProtID'].extend(prot_list)
            #  Records of sps: O.brachyantha
            if len(records) == 7:
                if records[6]:
                    gene_hash[gene_id]['ProtID'][records[6]] = '-'
#                    prot_list.append(records[6])
#                    gene_hash[gene_id]['ProtID'].extend(prot_list)
#     print(gene_hash)
    return gene_hash

    file_reader.close()

def geneParserPandas(infile):
    array = pd.read_csv(infile, sep="\t", delimiter=None, dtype='str')
    print array
    return array

def qtlParser(infile):
    headers = ['QTLid', 'Name', 'Symbol', 'TOid', 'Category', 'TraitName','TraitSymbol', 'Chromosome', 'Start', 'End']
    qtl_ds = list()
#    pp = pprint.PrettyPrinter(indent=4)
     
    fileHandle = open(infile, "r")
    lines = fileHandle.readlines()
    lines.pop(0) # remove header
    for line in lines:
        line = re.sub('\n$', '', line)
        items = line.split('\t')        
        qtl_ds.append(dict(zip(headers, items)))
    
    fileHandle.close()  
    return qtl_ds    
    
'''
def riceCyc(input_files):
    headers = ['Gene', 'Name', 'ReationId', 'ReactionName', 'EC', 'PathwayId', 'PathwayName']
    
    pw_ds = list()

    
    
    for input_file in input_files:
        fileHandle = open(input_file, "r")
    
#    print "*****Parsing RiceCyc data **********\n"
           
        lines = fileHandle.readlines()
        lines.pop(0)    
        for line in lines:
            line = re.sub('\n$', '', line)
            records = line.split('\t')
            pw_ds.append(dict(zip(headers, records)))
    
#    pw_ds.sort(key=lambda x: (x['PathwayId'], x['ReationId']))
    
    fileHandle.close()
    return pw_ds
    print "AraCyc data has been parsed!\n"
    print "**********************************\n\n"
'''

def CycParser(in_files):
    pp = pprint.PrettyPrinter(indent=4)
    pw_datastucture = {}
        
    for in_file in in_files:
        fileHandle = open(in_file, "r")
        lines = fileHandle.readlines()
        lines.pop(0)

        for line in lines:
#            ec_codes = {} #list()
#            reactions = {}
            pathways = {} 
            line = re.sub('\n$', '', line)
            records = line.split('\t')
            gene_id = records.pop(0)
            

#            reactions[records[1]] = {'ReactionName': records[2], 'EC': {records[3]: '-'}}
            
#            pathways[records[4]] = {'PathwayName': records[5], 'Reactions': reactions}
#            if records[4] not in pathways:
#                pathways[records[4]] = {'PathwayName': records[5], 'Reactions': {}}
#            if records[4] in pathways:
#                pathways[records[4]]['Reactions'].update(reactions)

            pathways[records[4]] = records[5]
            if gene_id in pw_datastucture:
                pw_datastucture[gene_id]['Pathways'].update(pathways)#= {records[4]: {'PathwayName': records[5], 'Reactions': []} }
                
            if gene_id not in pw_datastucture:
                pw_datastucture[gene_id] = {
                                            'Name': records[0],
                                            'Pathways': {}
                                            }
            if gene_id in pw_datastucture:
                pw_datastucture[gene_id]['Pathways'].update(pathways)#= {records[4]: {'PathwayName': records[5], 'Reactions': []} }
#            pp.pprint(reactions) #reactions pathways
    fileHandle.close()
    return pw_datastucture             


def getStrandValue(strandVar):
    if strandVar == "-1":
        strand = "-1"
        positionVar = "ReverseStrandPosition"
    else:
        strand = "1"
        positionVar = "ForwardStrandPosition"
    return (strand,positionVar)

''' 
 RDF Converters 
'''             
def grameneGeneRDF(files, output_dir,type='run'): #def grameneGeneRDF(files, output_dir): for test > type='test'
    rdf_buffer = ''
    output_file_name = os.path.split(os.path.splitext((files)[0])[0])[1]
    gene_counter = 0
    current_taxon_id = ''
    rdf_buffer = ''
    #genome_assembly = ''
#    pp = pprint.PrettyPrinter(indent=4)
    turtle_file = output_file_name + ".ttl"
    output_file = os.path.join(output_dir, turtle_file)
    output_opener = open(output_file, "w")
    chromosome_size = {'39947': {'size' : ['1:1-43270923:1', '2:1-35937250:1', '3:1-36413819:1', '4:1-35502694:1',
                                         '5:1-29958434:1', '6:1-31248787:1', '7:1-29697621:1', '8:1-28443022:1',
                                         '9:1-23012720:1', '10:1-23207287:1', '11:1-29021106:1', '12:1-27531856:1', 'Mt:1-402710:1', 'Pt:1-134481:1'],
                                 'genome_assembly' : 'IRGSP-1.0'},
                       '65489': {'size' : ['1:1-36915442:1', '2:1-31686972:1', '3:1-33311619:1', '4:1-27574323:1',
                                         '5:1-24206129:1', '6:1-25711811:1', '7:1-24128185:1', '8:1-22678724:1',
                                         '9:1-19219615:1', '10:1-19274048:1', '11:1-23014695:1', '12:1-20550741:1'],
                                 'genome_assembly': 'O.barthii_v1'},
                       '4533': {'size' : ['1:1-33916305:1', '2:1-27085147:1', '3:1-29620734:1', '4:1-21479432:1',
                                          '5:1-20131956:1', '6:1-21794218:1', '7:1-18603524:1', '8:1-18224760:1',
                                          '9:1-14107269:1', '10:1-14643570:1', '11:1-16001410:1', '12:1-15318893:1'],
                                 'genome_assembly' : 'Oryza_brachyantha.v1.4b'},
                       '4538':  {'size' : ['1:1-32613412:1', '2:1-29142869:1', '3:1-33053699:1', '4:1-26299011:1',
                                           '5:1-23192814:1', '6:1-24174485:1', '7:1-21799424:1', '8:1-20292731:1',
                                           '9:1-17607432:1', '10:1-16910673:1', '11:1-20796451:1', '12:1-19154523:1'],
                                 'genome_assembly' : 'Oryza_glaberrima_V1'},
                       '40149':  {'size' :[],'genome_assembly' : 'Oryza_meridionalis_v1.3'},  #  Oryza_meridionalis_v1.3
                       '4572' :  {'size' :[],'genome_assembly' : 'ASM34745v1'},  # t. urartu
                       '39946':  {'size' : ['1:1-47283185:1', '2:1-38103930:1', '3:1-41884883:1', '4:1-34718618:1',
                                            '5:1-31240961:1', '6:1-32913967:1', '7:1-27957088:1', '8:1-30396518:1',
                                            '9:1-21757032:1', '10:1-22204031:1', '11:1-23035369:1', '12:1-23049917:1'],
                                  'genome_assembly' : 'ASM465v1'},  # O.s indica
                       '40148': {'size': ['1:1-46529941:1', '2:1-38039368:1', '3:1-37594838:1', '4:1-32884216:1',
                                          '5:1-30429317:1', '6:1-31548187:1', '7:1-28350858:1', '8:1-27185770:1',
                                          '9:1-23490084:1', '10:1-22967634:1', '11:1-27098305:1', '12:1-26741765:1'],
                                 'genome_assembly': 'Oryza_glumaepatula_v1.5'},  # # Oryza_glumaepatula_v1.5
                       '4536': {'size': ['1:1-42845077:1', '2:1-35065507:1', '3:1-36134596:1', '4:1-28646061:1',
                                          '5:1-28014461:1', '6:1-29411429:1', '7:1-24717764:1', '8:1-26703665:1',
                                          '9:1-20407407:1', '10:1-21549876:1', '11:1-24378308:1', '12:1-20076173:1'],
                                 'genome_assembly': 'Oryza_nivara_v1.0'},  #
                       '4537': {'size': ['1:1-46096743:1', '2:1-39559433:1', '3:1-38925377:1', '4:1-33711903:1',
                                         '5:1-31082981:1', '6:1-34615992:1', '7:1-31244610:1', '8:1-29853973:1',
                                         '9:1-26294017:1', '10:1-25992119:1', '11:1-28494620:1', '12:1-27944835:1'],
                                'genome_assembly': 'Oryza_punctata_v1.2'},  #
                       '4558':  {'size' : ['1:1-80884392:1', '2:1-77742459:1', '3:1-74386277:1', '4:1-68658214:1',
                                           '5:1-71854669:1', '6:1-61277060:1', '7:1-65505356:1', '8:1-62686529:1',
                                           '9:1-59416394:1', '10:1-61233695:1'],  'genome_assembly' : 'Sorghum_bicolor_NCBIv3'},  # Sorghum_bicolor_NCBIv3
                       '4565':  {'size' : ['1A:1-594102056:1', '1B:1-689851870:1',  '1D:1-495453186:1', '2A:1-780798557:1',
                                           '2B:1-801256715:1', '2D:1-651852609:1' , '3A:1-750843639:1', '3B:1-830829764:1',
                                           '3D:1-615552423:1', '4A:1-744588157:1', '4B:1-673617499:1', '4D:1-509857067:1',
                                           '5A:1-709773743:1', '5B:1-713149757:1', '5D:1-566080677:1', '6A:1-618079260:1',
                                           '6B:1-720988478:1', '6D:1-473592718:1', '7A:1-736706236:1', '7B:1-750620385:1',
                                           '7D:1-638686055:1', 'Un:1-480980714:1'],   'genome_assembly' : 'IWGSC'},  # Triticum aevestivum
                       '4529':  {'size' : ['1:1-39866532:1', '2:1-33962743:1', '3:1-34446443:1', '4:1-30521686:1',
                                           '5:1-26652778:1', '6:1-27870862:1', '7:1-26200591:1' , '8:1-25958679:1',
                                           '9:1-20482102:1', '10:1-20731201:1', '11:1-27785585:1', '12:1-23561512:1'],
                                 'genome_assembly' : 'OR_W1943'},  # Oryza_rufipogon
                       '4577':  {'size' : ['1:1-307041717:1' , '2:1-244442276:1', '3:1-235667834:1', '4:1-246994605:1',
                                           '5:1-223902240:1', '6:1-174033170:1' , '7:1-182381542:1', '8:1-181122637:1',
                                           '9:1-159769782:1',  '10:1-150982314:1', 'Pt:1-140384:1', 'Mt:1-569630:1'],
                                 'genome_assembly' : 'B73_RefGen_v4'},  # Zea Mays
                       '3702':  {'size' : ['1:1-30427671:1', '2:1-19698289:1', '3:1-23459830:1', '4:1-18585056:1',
                                           '5:1-26975502:1'],  'genome_assembly' : 'TAIR10'},  # Arabidopsis
                       '4528': {'size': [], 'genome_assembly': 'O_longistaminata_v1.0'},  # Arabidopsis

                       }
    # Printing Prefixes
    output_opener.write(base + "\t" + "<" + base_uri + "> .\n")
    output_opener.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    output_opener.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    output_opener.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    output_opener.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    output_opener.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    output_opener.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n")
    output_opener.write(pr + "\t" + rapdb_gene_ns + "<" + rapdb_gene_uri + "> .\n")
    output_opener.write(pr + "\t" + msu_ns + "<" + msu_uri + "> .\n")
    output_opener.write(pr + "\t" + tair_l_ns + "<" + tair_l_uri + "> .\n")                
    output_opener.write(pr + "\t" + up_ns + "<" + uniprot + "> .\n\n")
    output_opener.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n\n")
    output_opener.write(pr + "\t" + ncbi_tax_ns + "<" + ncbi_tax_uri + "> .\n\n")
    output_opener.write(pr + "\t" + dc_ns + "<" + dc_uri + "> .\n\n")
    output_opener.write(pr + "\t" + faldo_ns + "<" + faldo_uri + "> .\n\n")
    output_opener.write(pr + "\t" + xsd_ns + "<" + xsd + "> .\n")
    output_opener.write(pr + "\t" + skos_ns + "<" + skos + "> .\n")
    output_opener.write(pr + "\t" + sio_ns + "<" + sio_uri + "> .\n")
    output_opener.write(pr + "\t" + chromosome_ns + "<" + chromosome_uri + "> .\n")
    output_opener.write(pr + "\t" + uniprot_ns + "<" + uniprot_uri + "> .\n")
    output_opener.write(pr + "\t" + ncbi_gene_ns + "<" + ncbi_gene_uri + "> .\n")
    output_opener.write(pr + "\t" + faldo_ns + "<" + faldo + "> .\n")
    output_opener.write(pr + "\t" + res_ns + "<" + resource + "> .\n")
    
    '''
    Ajout du prefix pour la release des donnees
    '''
    #output_opener.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")
    chromosome_number = 1
    for tax_id in taxon_ids:
        if output_file_name == taxon_ids[tax_id] or re.sub('_', ' ', output_file_name) == taxon_ids[tax_id]:
            current_taxon_id = tax_id
    genome_assembly = chromosome_size[current_taxon_id]['genome_assembly']
    if chromosome_size[current_taxon_id]['size']:
        print('test chromo')
        for ch_size in chromosome_size[current_taxon_id]['size']:
        # rdf_buffer += chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + str(
        #     chromosome_number) + ":1-" + str(
        #     ch_size) + ":1" + "\n"
            rdf_buffer += chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + ch_size + "\n"
            rdf_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + ncbi_tax_ns + current_taxon_id + " ;\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Chromosome" + " ;\n"
            rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + output_file_name + ":" + genome_assembly + ":" + ch_size + "\"@en ;\n"
            rdf_buffer += "\t" + dc_ns + "identifier " + "\t" + " \"" + current_taxon_id + ":" + genome_assembly + ":" + ch_size + "\" ;\n"
            rdf_buffer += "\t" + base_vocab_ns + "genomeAssembly " + "\t" + " \"" + genome_assembly + "\" .\n\n"
            chromosome_number += 1
    output_opener.write(rdf_buffer)

    for gene_file in files:

        rdf_buffer = ''
#        tigr_prefix = ''
#        rapdb_prefix = ''
#        geneId_prefix = ''
#        tair_prefix = ''

#        turtle_file = "gramene_" + output_file_name + "_genes" + ".ttl"
#        output_file = os.path.join(output_dir, turtle_file)
#        output_file = "gramene_genes.ttl"
#        output_opener = open(output_file, "w")
        
        print "*************Parsing %s genome data ***********\n" % (output_file_name)
        
        gene_ds = geneParser(gene_file)
        #gene_ds = geneParserPandas(gene_file)
#        pp.pprint(gene_ds)

#        print "%s data has been parsed!\n" % (output_file_name)
#        print "*************************************\n\n"
        #ict(gene_ds)
        print(str(len(gene_ds.keys())) + " Genes parsed")
        print "************* %s RDF conversion begins***********\n" % (output_file_name)
        reference_ch = ''
        for gene_id in gene_ds:
            rdf_buffer =''
            regex_ch =''
            chromosome_nb = 1
            gene_counter += 1
            (strand, position) = getStrandValue(gene_ds[gene_id]['Strand'])
            rdf_buffer += ensembl_ns + gene_id + "\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"

            for tax_id in taxon_ids:
                if output_file_name == taxon_ids[tax_id] or re.sub('_', ' ', output_file_name)== taxon_ids[tax_id]:
                    current_taxon_id = tax_id
                    rdf_buffer += "\t" + base_vocab_ns + "taxon" + "\t\t" + ncbi_tax_ns + current_taxon_id + " ;\n"

            for record_item in gene_ds[gene_id]:
                if record_item == 'Name':
                    if gene_ds[gene_id][record_item]:
                        rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (gene_ds[gene_id][record_item]) + " ;\n"
                if record_item == 'Description':
                    if gene_ds[gene_id][record_item]:
                        rdf_buffer += "\t" + dc_ns + "description" + "\t" + '"%s"' % (gene_ds[gene_id][record_item].replace("'", "")) + " ;\n"
                # if record_item == 'Chromosome':
                #     rdf_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t" + '"%s"' % (gene_ds[gene_id][record_item]) + " ;\n"
                # if record_item == 'Start':
                #     rdf_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + '"%s"' % (gene_ds[gene_id][record_item]) + " ;\n"
                # if record_item == 'End':
                #     rdf_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + '"%s"' % (gene_ds[gene_id][record_item]) + " ;\n"
                if record_item == 'Biotype':
                    rdf_buffer += "\t" + base_vocab_ns + "has_biotype" + "\t" + '"%s"' % (gene_ds[gene_id][record_item]) + " ;\n"
                if record_item == 'ProtID':
                    if gene_ds[gene_id][record_item]:
                        proteins = gene_ds[gene_id][record_item].keys()
                        for protein in proteins:
                            rdf_buffer += "\t" + base_vocab_ns + "encodes" + "\t" + up_ns + protein + " ;\n"
                if record_item == 'RapID':
                    if gene_ds[gene_id][record_item]:
                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + rapdb_gene_ns + gene_id + " ;\n"
                        rdf_buffer += "\t" + rdfs_ns + "seeAlso" + "\t" + rapdb_gene_ns + gene_id + " ;\n"
                if record_item == 'TairLocus':
                    if gene_ds[gene_id][record_item]:
                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + tair_l_ns + gene_ds[gene_id][record_item] + " ;\n"
                if record_item == 'TIGRlocus':
                    if gene_ds[gene_id][record_item]:
#                        tigr_prefix = pr + "\t" + tigr_ns + "<" + tigr_uri + "> .\n"
                        tigr_loci = gene_ds[gene_id][record_item].keys()
                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + msu_ns + tigr_loci[0].split('.')[0] + " ;\n"
                        # for locus in tigr_loci:
                        #     rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" + "\t" + msu_ns + locus + " ;\n"
                if record_item == 'Ontology':
                    # if gene_ds[gene_id][record_item]:
                    #     print('Ok ontology ...')
                        ont_terms = gene_ds[gene_id][record_item].items()
                        for ont in ont_terms:
                            ont_id = ont[1].replace(":", "_")
                            # print(ont[0])
                            if ont[0] == 'molecular_function':
                                rdf_buffer += "\t" + base_vocab_ns + "has_function" + "\t" + obo_ns + ont_id + " ;\n"
                                rdf_buffer += "\t" + base_vocab_ns + "go_term" + "\t" + obo_ns + ont_id + " ;\n"
                                # print('1 ...')
                            if ont[0] == 'biological_process':
                                rdf_buffer += "\t" + base_vocab_ns + "participates_in" + "\t" + obo_ns + ont_id + " ;\n"
                                rdf_buffer += "\t" + base_vocab_ns + "go_term" + "\t" + obo_ns + ont_id + " ;\n"
                                # print('2 ...')
                            if ont[0] == 'cellular_component':
                                rdf_buffer += "\t" + base_vocab_ns + "located_in" + "\t" + obo_ns + ont_id + " ;\n"
                                rdf_buffer += "\t" + base_vocab_ns + "go_term" + "\t" + obo_ns + ont_id + " ;\n"
                                # print('3 ...')
                            if ont[0] == 'plant_anatomy':
                                rdf_buffer += "\t" + base_vocab_ns + "expressed_in" + "\t" + obo_ns + ont_id + " ;\n"
                            if ont[0] == 'plant_structure_development_stage':
                                rdf_buffer += "\t" + base_vocab_ns + "expressed_at" + "\t" + obo_ns + ont_id + " ;\n"
                if record_item == 'StringID':
                    if gene_ds[gene_id][record_item]:
                        string_id = gene_ds[gene_id][record_item]
                        rdf_buffer += "\t" + base_vocab_ns + "has_dbxref" "\t" + "string:" + string_id + " ;\n"
                if record_item == 'Chromosome':
                    if gene_ds[gene_id][record_item]:
                        if chromosome_size[current_taxon_id]['size']:
                            for item in chromosome_size[current_taxon_id]['size']:
                                if gene_ds[gene_id]['Chromosome'] == item.split(':')[0]:
                                    reference_ch = item
                    # if not gene_ds[gene_id]['Chromosome']=='Mt':
                            chromosome_nb = gene_ds[gene_id]['Chromosome']
                            rdf_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  str(chromosome_nb) + ':' + str(gene_ds[gene_id]['Start']) + '-' + str(gene_ds[gene_id]['End']) + ":" + gene_ds[gene_id]['Strand'] +  " ;\n"
                        else:
                            reference_ch = str(chromosome_nb) + ':1-$:1'
                            rdf_buffer += "\t" + faldo_ns + "location" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                          str(chromosome_nb) + ':' + str(gene_ds[gene_id]['Start']) + '-' + str(
                                gene_ds[gene_id]['End']) + ":" + gene_ds[gene_id]['Strand'] + " ;\n"

            rdf_buffer = re.sub(' ;$', ' .\n\n', rdf_buffer)
            # Region
            rdf_buffer += chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['Start']) + '-' + str(gene_ds[gene_id]['End']) + ":" + gene_ds[gene_id]['Strand'] + "  \n"
            rdf_buffer += "\t" + rdfs_ns + "label" + "\t" + " \"" + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['Start']) + '-' + str(gene_ds[gene_id]['End']) + ":" + strand + "\";\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t" + faldo_ns + "Region" + " ;\n"
            rdf_buffer += "\t" + faldo_ns + "begin" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['Start'])  + ":" + strand + "  ;\n"
            rdf_buffer += "\t" + faldo_ns + "end" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['End']) + ":" + strand + "  .\n\n"

            # Position
            rdf_buffer += chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['Start'])  + ":" + strand
            rdf_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
            rdf_buffer += "  ;\n"
            rdf_buffer += "\t" + faldo_ns + "position" + "\t" + str(gene_ds[gene_id]['Start']) + " ;\n"
            rdf_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                              reference_ch + " .\n\n"

            # # Position 2
            rdf_buffer += chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                                  gene_ds[gene_id]['Chromosome']+ ':' + str(gene_ds[gene_id]['End']) + ":" + strand
            rdf_buffer += "\n" + "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + "ExactPosition" + " ;\n"
            rdf_buffer += "\t" + rdf_ns + "type" + "\t\t" + faldo_ns + position
            rdf_buffer += "  ;\n"
            rdf_buffer += "\t" + faldo_ns + "position" + "\t" + str(gene_ds[gene_id]['End']) + " ;\n"
            rdf_buffer += "\t" + faldo_ns + "reference" + "\t" + chromosome_ns + current_taxon_id + ":" + genome_assembly + ":" + \
                             reference_ch + " .\n\n"

            rdf_buffer = re.sub(' ;$', ' .\n\n', rdf_buffer)
            output_opener.write(rdf_buffer)
            if type == 'test': break

    output_opener.close()
#        if geneId_prefix:
#            output_opener.write(geneId_prefix)
#        if tigr_prefix:
#            output_opener.write(tigr_prefix)
#        if rapdb_prefix:
#            output_opener.write(rapdb_prefix)
#        if tair_prefix:
#            output_opener.write(tair_prefix)                            
        
    print "Number of genes in %s are: %s\n" % (output_file_name, str(gene_counter))
#    print "************* %s RDF completed ************\n" % (output_file_name)
#    print "Gramene gene data has been converted to RDF!\n"


def grameneQTLRDF(infile, output_dir):
    qtl_buffer = '' 
    to_hash = dict()
    qtl_counter = 0   
    turtle_file_name = "gramene.qtl.ttl"
    outfile = os.path.join(output_dir, turtle_file_name)
    outHandle = open(outfile, "w")
    
    print "*********** Parsing Gramene QTL data ***************\n"
    
    qtl_ds = qtlParser(infile)
    
#    print "Gramene QTL data has been parsed!\n"
#    print "*************************************\n" 
    
    print "************* Gramene QTL RDF conversion begins***********\n"
    
    outHandle.write(base + "\t" + "<" + base_uri + "> .\n")
    outHandle.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    outHandle.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    outHandle.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    outHandle.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    outHandle.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    outHandle.write(pr + "\t" + gr_qtl_ns + "<" + gramene_qtl + "> .\n\n")

    '''
    Ajout du prefix pour la realese des donnees
    '''
    #outHandle.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")


    for records in qtl_ds:
        qtl_buffer = ''
        qtl_counter += 1
        chrm = records['Chromosome'].replace("Chr. ", "")
        to_id = records['TOid'].replace(":", "_")
        
        qtl_buffer += gr_qtl_ns + records['QTLid'] + "\n"
        qtl_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "QTL" + " ;\n"
        #qtl_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
        #qtl_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + qtl_term + " ;\n"
        qtl_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (records['Name']) + " ;\n"
        if records['Symbol']:
            qtl_buffer += "\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['Symbol']) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "is_located_on" + "\t" + '"%s"' % (chrm) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "has_start_position" + "\t" + '"%s"' % (records['Start']) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "has_end_position" + "\t" + '"%s"' % (records['End']) + " ;\n"
        qtl_buffer += "\t" + base_vocab_ns + "has_trait" + "\t" + obo_ns + to_id + " ;\n"
            
#        if to_id not in to_hash:
#            outHandle.write(obo_ns + to_id + "\n")
#            outHandle.write("\t" + rdf_ns + "type" + "\t" + obo_ns + plant_trait_term + " ;\n") #base_vocab_ns + "Concept"
#            outHandle.write("\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + plant_trait_term + " ;\n")
#            outHandle.write("\t" + rdfs_ns + "label" + "\t" + '"%s"' % (records['TraitName']) + " ;\n")
#            outHandle.write("\t" + base_vocab_ns + "has_symbol" + "\t" + '"%s"' % (records['TraitSymbol']) + " ;\n")
#            outHandle.write("\t" + base_vocab_ns + "has_category" + "\t" + '"%s"' % (records['Category']) + " .\n")
#            to_hash[to_id] = 1
        qtl_buffer = re.sub(' ;$', ' .', qtl_buffer)
        outHandle.write(qtl_buffer)
    outHandle.close()
    print "Number of QTLs: %s\n" % (str(qtl_counter))
    print "********* Gramene QTL RDF completed ***********\n"

def CycRDF(data_stuc, output_dir):
    pw_hash = {}
#    react_hash = {}
#    gene_hash = {} 
#    previous_gene_id = ''
#    pw_buffer = ''
    gene_buffer = ''
#    react_buffer = ''
    pw_counter = 0
    gene_counter = 0
    tigr_pattern = re.compile(r'^LOC\_OS\d{1,2}G\d{5}\.\d$')
    sorghum_pattern = re.compile(r'^SB\d{2}G\d{6}\.\d$')
    alt_sorghum_match = re.compile(r'^SB\d{4}S\d{6}\.\d$')
    arabidopsis_pattern = re.compile(r'^AT[1-5]G\d{5}$')
    maize_pattern = re.compile(r'^GRMZM\d{1}G\d{6}')
    alt_maize_match = re.compile(r'^\w+\d{6}\.\d{1}\_\w+\d{3}')
    
    cyc_turtle = "gramene.cyc.ttl"
    outfile = os.path.join(output_dir, cyc_turtle)
    outputWriter = open(outfile, "w")
    
    print "*************Cyc RDF conversion begins***********\n"
    
    outputWriter.write(base + "\t" + "<" + base_uri + "> .\n")
    outputWriter.write(pr + "\t" + rdf_ns + "<" + rdf + "> .\n")
    outputWriter.write(pr + "\t" + rdfs_ns + "<" + rdfs + "> .\n")
    outputWriter.write(pr + "\t" + owl_ns + "<" + owl + "> .\n")
    outputWriter.write(pr + "\t" + base_vocab_ns + "<" + base_vocab_uri + "> .\n")
    outputWriter.write(pr + "\t" + obo_ns + "<" + obo_uri + "> .\n")
    outputWriter.write(pr + "\t" + swo_ns + "<" + swo_uri + "> .\n")
    outputWriter.write(pr + "\t" + sio_ns + "<" + sio_uri + "> .\n")
    outputWriter.write(pr + "\t" + pathway_ns + "<" + pathway_uri + "> .\n")
#    outputWriter.write(pr + "\t" + reaction_ns + "<" + reaction_uri + "> .\n")
    outputWriter.write(pr + "\t" + tigr_ns + "<" + tigr_uri + "> .\n")
    outputWriter.write(pr + "\t" + tigr_g_ns + "<" + tigr_g_uri + "> .\n")
    outputWriter.write(pr + "\t" + ensembl_ns + "<" + ensembl_plant + "> .\n\n")
#    outputWriter.write(pr + "\t" + ec_code_ns + "<" + ec_code_uri + "> .\n\n")

    '''
    Ajout du prefix pour la realese des donnees
    '''
    #outputWriter.write(pr + "\t" + res_ns + "<" + resource + "> .\n\n")


    # Genes
    for gene in data_stuc:
        gene_buffer = ''
        gene_counter += 1
        
        # Data from RiceCyc
        if tigr_pattern.match(gene):
            r_locus = list(gene)
            r_locus[5] = "s"
            r_locus[8] = "g"
            gene_locus = "".join(r_locus)            
            gene_buffer += tigr_ns + gene_locus + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + mrna_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "39947" + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "develops_from" + "\t" + tigr_g_ns + gene_locus[:14] + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]

        # Data from SorghumCyc        
        if sorghum_pattern.match(gene):
            s_g_id = list(gene)
            s_g_id[1] = "b"
            s_g_id[4] = "g"
            sorghum_gene_locus = "".join(s_g_id)            
            gene_buffer += ensembl_ns + sorghum_gene_locus + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + mrna_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "4558" + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "develops_from" + "\t" + ensembl_ns + sorghum_gene_locus[:11] + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]
        if alt_sorghum_match.match(gene):
            alt_s_g_id = list(gene)
            alt_s_g_id[1] = "b"
            alt_s_g_id[6] = "s"
            alt_sorghum_gene_locus = "".join(alt_s_g_id)
            gene_buffer += ensembl_ns + alt_sorghum_gene_locus + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + mrna_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "4558" + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "develops_from" + "\t" + ensembl_ns + alt_sorghum_gene_locus[:13] + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]

        #Data from MaizeCyc                    
        if maize_pattern.match(gene):
            gene_buffer += ensembl_ns + gene + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "4577" + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]
        if alt_maize_match.match(gene):
            gene_buffer += ensembl_ns + gene + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "4577" + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]

        # Data from AraCyc (Gramene)
        if arabidopsis_pattern.match(gene): 
            gene_buffer += ensembl_ns + gene + "\n"
            gene_buffer += "\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Gene" + " ;\n"
            #gene_buffer += "\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n"
            #gene_buffer += "\t" + rdfs_ns + "subClassOf" + "\t" + obo_ns + gene_term + " ;\n"
            gene_buffer += "\t" + rdfs_ns + "label" + "\t" + '"%s"' % (data_stuc[gene]['Name']) + " ;\n"
            gene_buffer += "\t" + base_vocab_ns + "taxon" + "\t" + obo_ns + "NCBITaxon_" + "3702" + " ;\n"
            for pw in data_stuc[gene]['Pathways']:
                gene_buffer += "\t" + base_vocab_ns + "is_agent_in" + "\t" + pathway_ns + pw + " ;\n"
                pw_hash[pw] = data_stuc[gene]['Pathways'][pw]
        
        gene_buffer = re.sub(' ;$', ' .', gene_buffer)        
        outputWriter.write(gene_buffer)

    #Pathways
    for pw_id in pw_hash:
        pw_counter += 1
        outputWriter.write(pathway_ns +  pw_id + "\n")
        outputWriter.write("\t" + rdf_ns + "type" + "\t" + base_vocab_ns + "Pathway_Identifier" + " ;\n")
        #outputWriter.write("\t" + rdf_ns + "type" + "\t" + owl_ns + "Class" + " ;\n")
        #outputWriter.write("\t" + rdfs_ns + "subClassOf" + "\t" + sio_ns + met_pw_sio_term + " ;\n")
        #outputWriter.write("\t" + rdfs_ns + "subClassOf" + "\t" + swo_ns + biocyc_pw_term + " ;\n")
        outputWriter.write("\t" + rdfs_ns + "label" + "\t" + '"%s"' % (pw_hash[pw_id]) + " .\n")
        
    outputWriter.close()
    
    print "Number pathways and genes: %s and %s\n" % (str(pw_counter), str(gene_counter))
    print "******* Cyc RDF completed **********\n" 

    
'''
 Tools                    
'''        
def removeDuplicates(in_list):
    newlist = list(set(in_list))
    return newlist

ROOT_DIR='/Users/plarmande/Downloads/data/'

# ROOT_DIR='/Volumes/LaCie/AGROLD/agroLD_data_update_mai_2017'
gramene_genes_files = [ROOT_DIR + 'Oryza_meridionalis.txt']
gramene_genes_out =  '/Users/plarmande/Downloads/data/'
# gramene_qtl_out = ROOT_DIR + '/rdf/gramene_qtl_ttl/'


pp = pprint.PrettyPrinter(indent=4)

gramene_genomes = gramene_genes_files #gramene_genes_files
# g_parse = grameneParsers#oryzaBaseParser
print "***************** Gramene Genes data ********************\n"
print gramene_genes_files
# g_parse = geneParser(gramene_genomes)
# input_f = '/Oryza brachyantha.txt' #Oryza_barthii.txt' Oryza_sativa_japonica
# input_f = '/home/venkatesan/workspace/explore/test_files/gramene_genes/Oryza_brachyantha.txt' #Oryza_barthii.txt' Oryza_sativa_japonica
#geneHash = geneParser(gramene_genomes)#grameneQTLRDF(gramene_qtl_dir, gramene_qtl_out) oryzaBaseRDF(oryzabase_file, oryzaBase_output) grameneGeneRDF(gramene_genomes, gramene_genes_out)
#pp.pprint(geneHash)
grameneGeneRDF(gramene_genomes, gramene_genes_out)
print "********************************************************\n\n"

