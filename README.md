
# AgroLD ETL

AgroLD is a RDF knowledge base that consists of data integrated from a variety of plant resources and ontologies. AgroLD ETL is the Python packages developed to transform plant datasets in RDF. Packages are developped for data standards such as GFF,GAF, VCF and specific plant databases.


# Contact

* pierre.larmande_at_ird.fr


# Valorization

* [https://doi.org/10.1371/journal.pone.0198270](https://doi.org/10.1371/journal.pone.0198270)
* [https://tel.archives-ouvertes.fr/IBC/hal-01176903v1](https://tel.archives-ouvertes.fr/IBC/hal-01176903v1) 
* [http://ceur-ws.org/Vol-1546/poster_55.pdf](http://ceur-ws.org/Vol-1546/poster_55.pdf)



# Contributing

* Intellectual property belongs to IRD, CIRAD, IBC, INRA, IFB, ELIXIR, and SouthGreen development platform.
* Written by Aravind Venkatesan, Gildas Tagny Ngompe, Nordine El Hassouni, Imene Chentli, Valentin Guignon, Clement Jonquet, Manuel Ruiz, Pierre Larmande. 
* Copyright 2014-2018


# The Architecture of project


AgroLD project is composed of two component: [AgroLD_ETL](/AgroLD_ETL)  and  [agrold_webapp](/agrold_webapp)


* The first component is a set of Parser and wrapper for translate a dataset. Follow this link for to know what data have been translated in RDF: [Documentation](http://volvestre.cirad.fr:8080/agrold/documentation.jsp)  

* The second component is the web application who is connected at the triple store for to make queries.
 
***

```
 AgroLD_ETL
	-> model
	-> rdf_ttl
	-> riceKB
	-> test_files
	-> riceKBpipeline.py
```
***


### AgroLD_ETL


Contains parsers and model used to convert data considered for AgroLD to RDF.

* [model](/AgroLD_ETL/model): All documents who describe how data are transformed
* [rdf_ttl](/AgroLD_ETL/rdf_ttl): All output of transformation sort by dataset
* [riceKB](/AgroLD_ETL/riceKB): Contains scripts used for each data set
* [test_files](/AgroLD_ETL/test_files): All test files in input ( heterogeneous format: csv, tabbed files, gff3 )
* [riceKBpipeline.py](/AgroLD_ETL/riceKBpipeline.py): Script file where we have centralised all execution


The type of each dataset is different, GFF, HapMap, CSV and VCF. In first time we have developed parser for build a dictonary, 
because is easy to browse a dictionary and create the RDF 


# How to use AgroLD_ETL

For example if you want to execute a gff3 parser for to build a dictionary.
define a input file

```
AgroLD/AgroLD_ETL/riceKB/gffParser.py

path = '/os_file_gff3/file.gff3'     # The input

ds = parseGFF3(path)   # The parsing file

```

> **The dictionary :** The GFF3 Parser is a generic fonction who build a dictionary, it easy to browse this dictionary for build RDF data 


```
{   'attributes': {   'Dbxref': 'InterPro:IPR005333,InterPro:IPR017887',
                          'ID': 'BGIOSGA000770-TA',
                          'Name': 'BGIOSGA000770-TA',
                          'Parent': 'BGIOSGA000770'},
        'end': 35414873,
        'phase': None,
        'score': None,
        'seqid': 'Osi01',
        'source': 'glean',
        'start': 35413950,
        'strand': '-',
        'type': 'mRNA'},
    {   'attributes': {   'Derives_from': 'BGIOSGA000770-TA',
                          'ID': 'BGIOSGA000770-TA_protein',
                          'Name': 'BGIOSGA000770-TA'},
        'end': 35414873,
        'phase': None,
        'score': None,
        'seqid': 'Osi01',
        'source': 'glean',
        'start': 35413950,
        'strand': '-',
        'type': 'polypeptide'},
    {   'attributes': {   'Parent': 'BGIOSGA000770-TA'},
        'end': 35414873,
        'phase': '0',
        'score': None,
        'seqid': 'Osi01',
        'source': 'glean',
        'start': 35413950,
        'strand': '-',
        'type': 'CDS'},
    {   'attributes': {   'Parent': 'BGIOSGA000770-TA'},
        'end': 35414873,
        'phase': None,
        'score': None,
        'seqid': 'Osi01',
        'source': 'glean',
        'start': 35413950,
        'strand': '-',
        'type': 'exon'},

```

### Documentation

- AgroLD includes data on the following species on :  [Species](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#species)
- Ontologies in AgroLD : [Ontologies](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#ontologies)
- Data sources in AgroLD : [Data](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#sources)
- Species specific break down of the data sources : [Link](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#break-down)
- Graph Names : [Link](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#graphs)
- URIs :  [Link](http://volvestre.cirad.fr:8080/agrold/documentation.jsp#uri)



