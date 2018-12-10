# GreenPhyl RDF data model

## Prefixes: 
 * @base <http://www.southgreen.fr/agrold/> . 
 * @prefix rdf:<http://www.w3.org/1999/02/22­rdf­syntax­ns#> . 
 * @prefix rdfs:<http://www.w3.org/2000/01/rdf­schema#> .
 * @prefix xsd:<http://www.w3.org/2001/XMLSchema#> . 
 * @prefix agrold_vocabularyulary:<http://www.southgreen.fr/agrold/vocabulary/> .
 * @prefix obo:<http://purl.obolibrary.org/obo/> .
 * @prefix uniprot:<http://purl.uniprot.org/uniprot/> .
 * @prefix tigrlocus:<h​ttp://www.southgreen.fr/agrold/tigr.locus/>​.
 * @prefix tairlocus:<h​ttp://identifiers.org/tair.locus/>​.
 * @prefix pubmed:<h​ttp://identifiers.org/pubmed/>​. 
 * @prefix interpro:<h​ttp://identifiers.org/interpro/>​. 
 * @prefix ncbi_taxon:<h​ttp://purl.obolibrary.org/obo/NCBITaxon_>​.
 * @prefix  greenphyl_family:<h​ttp://www.southgreen.fr/agrold/greenphyl.family/>​.
 * @prefix  greenphy_annotation:<h​ttp://www.southgreen.fr/agrold/greenphyl.annotation/>​.

## Predicates to be used:

* rdf:type
* rdf:subject
* rdf:predicate
* rdf:object
* rdfs:subClassOf
* rdfs:label
* agrold_vocabulary:description
* agrold_vocabulary:taxon
* agrold_vocabulary:contains
* agrold_vocabulary:is_orthologous_to 
* agrold_vocabulary:is_paralogous_to 
* agrold_vocabulary:is_member_of
* agrold_vocabulary:xRef
* agrold_vocabulary:has_annotation
* agrold_vocabulary:has_score
* agrold_vocabulary:assigned_by 
* agrold_vocabulary:has_uniprot_accession 
* agrold_vocabulary:number_of_sequences 
* agrold_vocabulary:has_thresold 
* agrold_vocabulary:curation_status 
* agrold_vocabulary:has_uniprot_accession 
* agrold_vocabulary:has_go_identifier
* obo:BFO_0000056 # participates_in ­ Biological process
* obo: BFO_0000085 # has_function ­ Molecular Function 
* obo:BFO_0000082 # located_in ­ Cellular Component


## Protein information:
* ID: Uniprot <h​ttp://purl.uniprot.org/uniprot/>​
* Name: Literal
* Description: literal
* Protien Ontology term: SO_0000104 <obo:SO_0000104>
* Taxon: NCBI <http://purl.obolibrary.org/obo/NCBITaxon_>
* Gene: MSU/TIGR Locus ID <h​ttp://www.southgreen.fr/agrold/tigr.locus/>​(O.sativa); 
* TAIR ID <h​ttp://identifiers.org/tair.locus/>​(A.thaliana)
* IPR: InterPro <h​ttp://identifiers.org/interpro/>​ 
* GO: ID <h​ttp://purl.obolibrary.org/obo/>​ 
* Homology info: 
* Orthology: Uniport ID (URI)
* Paralogy: Uniport ID (URI)
* Score: literal 
* Annotation: Literal
* DB xref: URI or literal, depends on the kind of xref e.g IDs 
* Reference: Pubmed ID <h​ttp://identifiers.org/pubmed/>​


## Protein RDF Model:

> E.g. Arabidopsis
> **Note:​**For the corresponding ortholog protein a prefix need to be based on the source organism.

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

* ID : Internal <h​ttp://www.southgreen.fr/agrold/greenphyl.family/>​ 
* cluster (family): OBI_0000251
* Name: Literal
* Level: Literal
* Number_of_sequences: Literal
* Evidence:Literal
* GO: ID <h​ttp://purl.obolibrary.org/obo/>​
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
	* agrold_vocabulary:curation_status “Status” ; # A​nnotation in progress 
	* agrold_vocabulary:xRef pubmed:XXXXXX .
  