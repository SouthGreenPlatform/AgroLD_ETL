# GreenPhyl RDF data model

## Prefixes: 
 * @base <http://www.southgreen.fr/agrold/> . 
 * @prefix rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>. 
 * @prefix rdfs:<http://www.w3.org/2000/01/rdf-schema#>.
 * @prefix xsd:<http://www.w3.org/2001/XMLSchema#>. 
 * @prefix agrold_vocabulary:<http://www.southgreen.fr/agrold/vocabulary/>.
 * @prefix obo:<http://purl.obolibrary.org/obo/>.
 * @prefix uniprot:<http://purl.uniprot.org/uniprot/>.
 * @prefix tairlocus:<http://identifiers.org/tair.locus/>.
 * @prefix pubmed:<http://identifiers.org/pubmed/>. 
 * @prefix interpro:<http://identifiers.org/interpro/>. 
 * @prefix ncbi_taxon:<http://purl.obolibrary.org/obo/NCBITaxon_>.
 * @prefix greenphyl_family:<http://www.southgreen.fr/agrold/greenphyl.family/>.
 * @prefix greenphy_annotation:<http://www.southgreen.fr/agrold/greenphyl.annotation/>.
 
### New prefixes :
* @prefix skos:<http://www.w3.org/2004/02/skos/core#>.
* @prefix owl:<http://www.w3.org/2002/07/owl#>.
* @prefix dcterms:<http://purl.org/dc/terms/>.
* @prefix dc:<http://purl.org/dc/elements/1.1/>.
* @prefix doi:<http://dx.doi.org/>.
* @prefix faldo:<http://biohackathon.org/resource/faldo#>.
* @prefix agrold_resource:<http://www.southgreen.fr/agrold/resource/>.
* @prefix tigr_gene:<http://identifiers.org/ricegap/>
* @prefix ensembl_plant:<http://identifiers.org/ensembl.plant/>.
* @prefix sio:<http://semanticscience.org/resource/> .
* @prefix taxon:<http://identifiers.org/taxonomy/>

### Deprecated prefixes:
* @prefix tigrlocus:<http://www.southgreen.fr/agrold/tigr.locus/>.


## Predicates to be used:

* rdf:type
* rdf:subject
* rdf:predicate
* rdf:object
* rdfs:subClassOf
* rdfs:label
* agrold_vocabulary:has_score
* agrold_vocabulary:assigned_by 
* agrold_vocabulary:has_uniprot_accession 
* agrold_vocabulary:number_of_sequences 
* agrold_vocabulary:has_thresold 
* agrold_vocabulary:curation_status
* obo:BFO_0000056 # participates_in - Biological process
* obo: BFO_0000085 # has_function - Molecular Function 
* obo:BFO_0000082 # located_in - Cellular Component

### New predicates: 
* obo:RO_0002162 # in taxon
* dcterms:identifier
* dc:description
* rdfs:comment
* rdfs:seeAlso
* owl:sameAs
* obo:RO_0002162 # in taxon
* obo:RO_0002350 # is member of
* obo:RO_0002524 # contains
* sio:SIO_000255 # has annotation
* sio:SIO_010082 # is translated into
* sio:SIO_010083 # is translated from 
* sio:SIO_010081 # is transcribed from
* sio:SIO_010080 # is transcribed into
* sio:SIO_010078 # encodes
* sio:SIO_000558 # is orthologous to
* sio:SIO_000630 # is paralogous to
* skos:altLabel
* skos:prefLabel
* skos:prefSymbol
* skos:altSymbol
* agrold_vocabulary:classifiedWith
* dcterms:references

## Deprecated predicates
* agrold_vocabulary:xRef
* agrold_vocabulary:description
* agrold_vocabulary:taxon
* agrold_vocabulary:is_orthologous_to 
* agrold_vocabulary:is_paralogous_to 
* agrold_vocabulary:is_member_of
* agrold_vocabulary:contains
* agrold_vocabulary:has_annotation
* agrold_vocabulary:has_go_identifier

## Protein information:
* ID: Uniprot <http://purl.uniprot.org/uniprot/>
* Name: Literal
* Description: literal
* Protien Ontology term: SO_0000104 <obo:SO_0000104>
* Taxon: NCBI 	<http://identifiers.org/taxonomy/>
* Gene: MSU/TIGR Locus ID <http://www.southgreen.fr/agrold/tigr.locus/>(O.sativa); 
* TAIR ID <http://identifiers.org/tair.locus/> (A.thaliana)
* IPR: InterPro <http://identifiers.org/interpro/> 
* GO: ID <http://purl.obolibrary.org/obo/> 
* Homology info: 
* Orthology: Uniport ID (URI)
* Paralogy: Uniport ID (URI)
* Score: literal 
* Annotation: Literal
* DB xref: URI or literal, depends on the kind of xref e.g IDs 
* Reference: Pubmed ID <http://identifiers.org/pubmed/>


## Protein RDF Model:

> E.g. Arabidopsis
> **Note: **For the corresponding ortholog protein a prefix need to be based on the source organism.

* tair:ATXXXX
	* rdf:type 							agrold_vocabulary:Protein ; 
	* rdfs:subClassOf 					obo:SO_0000104 ;
	* rdfs:label 							“name” ; 
	* agrold_vocabulary:description 		“description” ; 
	* agrold_vocabulary:taxon 			ncbi_taxon:XXXX ; 
	* obo:BFO_0000056 # participates_in	obo:GO_ID ;
	* obo: BFO_0000085 # has_function		obo:GO_ID ; 
	* obo:BFO_0000082 # located_in		obo:GO_ID ; 
	* agrold_vocabulary:contains 			interpro:XXXXX ;
	* agrold_vocabulary:is_member_of 		greenphyl_family:XXXXX ; # family ID
	* agrold_vocabulary:is_orthologous_to	uri:XXXXXX ;
	* agrold_vocabulary:has_annotation  	greenphy_annotation:ID1XXXXX_ID2XXXX ;
	* agrold_vocabulary:is_paralogous_to	uri:XXXXXX ;
	* agrold_vocabulary:has_annotation	greenphy_annotation:ID1XXXXX_ID3XXXX ;
	* agrold_vocabulary:has_uniprot_accession uniprot:XXXXX ; 
	* agrold_vocabulary:xRef 				pubmed:XXXXXX ;
	* agrold_vocabulary:xRef 				tairlocus:XXXXXX .


## Reification 

* greenphy_annotation:ID1XXXXX_ID2XXXX
	* rdf:type	rdf:Statement ;
	* rdf:subject	tair:XXXXX ; # ID1
	* rdf:Predicate agrold_vocabulary:is_orthologous_to ; 
	* rdf:object tair:XXXXX ; # ID2
	* agrold_vocabulary:has_score “score” ; # Literal 
	* agrold_vocabulary:assigned_by “GreenPhyl” .




## Family information:

* ID : Internal <http://www.southgreen.fr/agrold/greenphyl.family/> 
* cluster (family): OBI_0000251
* Name: Literal
* Level: Literal
* Number_of_sequences: Literal
* Evidence:Literal
* GO: ID <http://purl.obolibrary.org/obo/>
* DB xref: URI or literal, depends on the kind of xref e.g IDs


## Family RDF model:

* greenphyl_family:GXXXX
	* rdf:type agrold_vocabulary:Protein_Family ;
	* rdf:subClassOf obo:OBI_0000251 ;
	* rdfs:label “name” ;
	* agrold_vocabulary:description “description” ; 
	* agrold_vocabulary:has_go_identifier obo:GO_ID ; 
	* agrold_vocabulary:evidence “Evidence” ; 
	* agrold_vocabulary:number_of_sequences “#”^^xsd:integer ; 
	* agrold_vocabulary:has_thresold “#”^^xsd:integer ; # e.g. 1, 2 ,3 
	* agrold_vocabulary:curation_status “Status” ; # Annotation in progress 
	* agrold_vocabulary:xRef pubmed:XXXXXX .
  
