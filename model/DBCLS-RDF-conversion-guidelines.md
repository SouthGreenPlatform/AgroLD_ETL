# DBCLS database RDF conversion guidelines (translated from [Japanese version](https://github.com/dbcls/rdfizing-db-guidelines/blob/master/dbcls-rdfizing-db-guidelines-ja.md) with Google)

## Introduction - Linked Data Concept

[RDF](http://www.w3.org/TR/rdf11-concepts/) is a framework that describes data with low ambiguity and high machine readability. However, guidelines are not provided as to how RDF can be used to actually describe the data to be described, whether the content desired to be expressed can be described, and whether it will be better from the viewpoint of subsequent use. This is a hurdle for beginners, especially when converting databases to RDF. Fortunately, [BioHackathon](http://biohackathon.org/) , [SPARQLthon](http://wiki.lifesciencedb.jp/mw/SPARQLthon), and so on have also had various discussions and knowledge accumulated. I think that it is time to come up with guidelines to guide guidelines when converting data into RDF in this community. The goal of this guideline is to make it possible to create an RDF that reduces the burden of RDF conversion work and can be integrated appropriately with other data by referring to it.

As a basic spirit, it is recommended to create data based on the [Linked Data concept](http://www.w3.org/DesignIssues/LinkedData.html) by Tim Berners-Lee . In the Linked Data concept,

1. Use URIs as names for things → Name and name items using URI
2. Use HTTP URIs so that people can look up those things. → Use URIs starting with http: // as names that can be accessed by widely popular software so that users can look over them (widely spread Because there are some URIs that can not be accessed by software you are using)
3. When someone looks up a URI, provide useful information, using the standards (RDF *, SPARQL) → Provide useful information according to standards such as RDF and SPARQL when accessing URI
4. Include links to other URIs. So that they can discover more things. → Include links to other URIs so that you can trace more information

The four principles are advocated.

RDF is an abbreviation for [Resource Description Framework](https://en.wikipedia.org/wiki/Resource_Description_Framework) , and URI is an abbreviation for [Uniform Resource Identifier](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier) . RDF is a framework (standard) for representing information on [resources](https://en.wikipedia.org/wiki/Web_resource) pointed to by a URI on the web, not a specific file format. In fact, there are several formats that describe information using RDF, such as RDF / XML, N-Triples, Turtle, JSON-LD, etc. Turtle format is highly readable not only from machines but also from people, Because there are many points in common with the notation method, in this guideline, we will make notation according to Turtle format.

### 1 Guidelines for making RDF of database contents

Below, from the findings accumulated so far, the recommended guidelines for RDF conversion of databases are listed.

#### 1.1 Guidelines for Designing URIs

##### 1.1.1 Use highly persistent URI

Guideline of using a highly persistent URI is, according to Tim Berners-Lee " [Cool URI does not change](https://www.w3.org/Provider/Style/URI) " has begun from the declaration that, [cool URI for from the W3C](https://www.w3.org/TR/cooluris/) is released.

RDF data describes information about URI resources. In Principle 2 of Linked Data, it is recommended to use URI starting with http: // and in Principle 3 information on URI resources can be obtained using standards such as RDF by accessing the URI We recommend. However, if the URI that was originally accessible changes to a different URI for some reason or it disappears, the value of RDF data significantly decreases. For example, many of life science databases are made by research institutes, but it is not unusual for URIs of organizations to change due to organization reorganization, integration, abolition etc. Therefore, if RDF data is described by URI using domain name of the organization, even if the database is to be operated with a new organization luckily, the URI can not be accessed from the RDF data. In order to avoid such a situation, it is recommended to consider the following two methods.

* Get your own domain

As a common case, there is a way to get your own domain like "database name .org" like uniprot.org. Thus, even if the entity operating the database changes, the same URI can be used continuously. The disadvantage is that there are other persistence problems such as costly to maintain your own domain and who maintains that domain.

* Using PURL

Another commonly used case is to describe data using [Persistent URL (PURL)](https://en.wikipedia.org/wiki/Persistent_uniform_resource_locator) which is a [persistent URL](https://en.wikipedia.org/wiki/Persistent_uniform_resource_locator) . RDF is described by using URI of PURL, and by transferring (redirecting) data from PURL to the URL of the entity operating data, even when the entity is transferred, the same URI is maintained without changing the RDF data You can use it. [find service providing PURL] By managing the transfer with the PURL service, you can always point to the latest URI that you can access at the moment.

##### 1.1.2 A URI indicating a resource describes an ID for identifying it in the end of the URI

For the path part of the URI following the highly persistent domain part, it is recommended that the ID that identifies the description subject resource representing mono or koto is at the end of the URI, and the slash (/) is used immediately before that It will be. For example UniProt assigns a URI of http://purl.uniprot.org/uniprot/Q6GZX3 for proteins identified by the ID Q6GZX3 . By configuring URI like this, it is technically realized to construct a service "Provide useful information according to standards such as RDF and SPARQL when accessing URI" indicated in Linked Data Principle 3 It becomes easy.


In the ontology URI, when creating individual pages for each ID indicating concepts that make up the ontology, a slash (/) is provided immediately before the ID and the entire ontology is provided in one page as described above. If you do, it is a good idea to have a hash (#) in front of the ID instead of a slash. As an example, since UniProt's ontology '''http://purl.uniprot.org/core/provides''' concepts such as Protein on individual pages, '''http://purl.uniprot.org/core/Proteinslash''' (/) is used as in UniProt ontology . On the other hand, in the ontology [FALDO](http://biohackathon.org/resource/faldo) of the array position information, the concept such as Position is provided on the same page, '''http://biohackathon.org/resource/faldo#Positionso''' a hash (#) is used as with the URI fragment notation 

##### 1.1.3 Handling of ID including version information

In NCBI Gene and NCBI Protein etc, ID including version number and ID not included are provided (eg [https://www.ncbi.nlm.nih.gov/protein/NP_003024.1](https://www.ncbi.nlm.nih.gov/protein/NP_003024.1) and [https://www.ncbi.nlm.nih.gov/protein/NP-003024](https://www.ncbi.nlm.nih.gov/protein/NP-003024)). When using such an ID, which ID should be used is a problem. When attaching importance to links with other data sets, linking between data sets will be easier to maintain if IDs that do not include version numbers are used. On the other hand, if you want to describe information about a specific version of the ID, you should use the ID that contains the version number. That is, you need to use a more appropriate ID depending on the situation. However, from the viewpoint of maintaining links between datasets, it is preferable to provide a link to an ID without a version number even if you use an ID with a version number. If you want to know the latest version, [TogoWScan](http://togows.org/) use the following to obtain the latest version of each database for the corresponding database in the following way.

'''http://togows.org/entry/ncbi-protein/145579718/version'''

#### 1.2 Guidelines for creating RDF
##### 1.2.1 URI resource is defined as an instance of class of ontology

'''rdf:type''' It is important that ontology classes are specified in to clearly and concisely express what the resource indicated by URI means . In particular, it is recommended to specify classes on the basis of ontology for the main resource.

'''uniprot:Q6GZX3 rdf:type core:Protein .'''

In the example above, [core:Protein](http://purl.uniprot.org/core/Protein) you can see that the Q6GZX3 entry of UniProt is . In addition to clarifying what the resource points, you can also rdf:typeperform efficient data search by specifying it also in the SPARQL search . Ideally it is possible to refer to the definition of ontology by accessing URI, but that is not essential. However, in natural languages, it is common for different contexts to have different meanings if different contexts are used. Describing the definition intended by the data creator in the reference destination of the URI is an effective means to avoid misunderstanding of the data user. This is a great advantage in describing scientific data.

The problem here is that, "Which class should I type in my resource?" If you find a suitable class for an existing ontology you should use it, but in many cases you will not find it. In that case, you need to create a simple ontology for your own data. UniProt and EBI RDF also define required classes and provide them with data.

In RDF,

* '''rdfs:Class''' When it is '''rdf:type''' the object of, that the subject is a class

	'''Myclass rdf:type rdfs:Class .'''

* When a URI of '''rdf:type''' a certain class is the object of the object, that the subject is an instance

	'''ex:111 rdf:type ex:Myclass .'''

Even if a subject is defined as a subclass of a class by, rdfs:Classeven if it is not explicitly typed '''rdfs:subClassOf''', its subject is also a class). However, to instantiate a URI as a class and as an instance complicates the semantic processing of data, please avoid it if it can be avoided.

'''#not recommended 
ex:111 rdf:type rdfs:Class # ← ex:111 はクラスでもあり  
ex:111 rdf:type ex:Myclass # ← インスタンスでもある、ことになる  '''

#####  1.2.2 Label URI Resources

In order for human beings to easily understand what the resource indicated by URI is, rdfs:labelit is useful that labels are written in more natural language. In particular, it is recommended to describe a concise label on the main resource. This is not related to enhancing machine readability, but it improves human readability, which makes it convenient for creating applications and for displaying the results of SPARQL search in an easy-to-read manner.

	'''uniprot:P51028 rdfs:label `“`wnt8a`”`@en .`'''

Also, as in this example, attaching a language tag is useful especially for multilingualization. If it is not English native, you may want to add a label even in your native language, or you may want to add comments. By attaching language tags, you can explicitly write in multiple languages, and you can easily do things like using only English labels when using in applications.

	'''mpo:MPO_03001 rdfs:label "Thermophilic"@en ,"高熱性"@ja .'''

I think that it is convenient to keep only one label in one language because it does not bother trouble when using it. If you want to give a different label, you can use rdfs:labelinstead '''skos:altLabel'''. Also, if there are multiple labels due to '''rdfs:label''', '''skos:altLabel''' etc., you '''skos:prefLabel''' can use one if you want to specify one of them as a representative label .


	
	'''chebi:CHEBI_17234 rdfs:label "D-Glucose"@en ;  
                  rdfs:label "D-グルコース""@ja ;  
                  skos:altLabel "Dextrose"@en ;  
                  skos:altLabel ""ブドウ糖""@ja ;  
                  skos:prefLabel "D-Glucose"@en .  '''

##### 1.2.3 Append ID label to URI resource

The URI itself functions as a global ID. However, URIs are symbolic and long as strings, so it is not suitable for humans to see when displaying SPARQL results. If the end of the URI contains a database-specific (local) ID, you can get the ID as a character string by cutting out the last / / of the URI. To do this, use SPARQL It takes extra time and effort to apply string processing with. Therefore, dcterms:identifierit is recommended to describe the ID string of the main resource with .

	'''uniprot:P51028 dcterms:identifier "P51028" .'''

Please '''dcterms:identifier''' do not declare multiple IDs using properties for a resource identified by ID.

	'''#The following are deprecated.  
	pdb:2RH1 dcterms:identifier "2RH1" .`  
	pdb:2RH1 dcterms:identifier "2rh1" .`'''

Although it is not possible to give a different ID from the original meaning of the identifier of ID, if you want to declare multiple identifiers for some reason (for example, in case of both upper case and lower case as in the example above, corresponded to inquiries Please declare using another property for case etc.).

	'''#There is no problem below.  
	pdb:2RH1 dcterms:identifier "2RH1" .  
	pdb:2RH1 skos:altLabel "2rh1" .   '''


##### 1.2.4 Add a link to another data set

In RDF, as shown in Principle 4 of Linked Data, by pasting a reference link to an external resource, the Web of data (Linked Data) is realized. In the case where the cross reference destination is a database entry, it is common to link to the URL that the reference DB entries can actually see. I call this a polite URL from the meaning of respecting the original site. However, as a problem of using the original URL,

* Even if the URL of the linked page is changed, the change will not be automatically applied to RDF once it is released
* There are cases where the original URL is not a cool URI (it is not suitable as an ID, such as a page dynamically generated by taking CGI arguments)

There are. In addition, there are cases where there are multiple URLs for referring to the same item (for example, [Taxonomy ID](http://info.identifiers.org/taxonomy/9606) etc is that sites such as NCBI, EBI, UniProt, Bio 2 RDF etc. provide different URLs, all of which are widely used 

Move
* Providing RDF using different URIs for individual databases to point to the same thing makes integrated searching impossible

There is also a problem.

As a method to avoid such problems,It is possible to use the URI of [http://identifiers.org/](http://identifiers.org/) . Identifiers.org is similar to PURL, but it is a service specialized for URI of the database, it can be selected if there are multiple destinations, and database metadata is provided by RDF etc. There are advantages. Identifiers.org is maintained by the life science community and can be added by request even if the database you want to see is not registered.

For this reason, when describing the cross reference,

	'''example)  
	ex:111 rdfs:seeAlso <http://pfam.xfam.org/family/PF01590> .  
	ex:111 rdfs:seeAlso <http://identifiers.org/pfam/PF01590> . '''

By attaching a reference link to both the polite URL and the URI of Identifiers.org, you can create an RDF that is easy to connect with external resources.

For NBDC RDF portal, it is recommended to link against the URI used in the RDF and the URI of Identifiers.org (if available) for the database for which RDF is published. As for the URI of Identfiers.org, we plan to generate a triple on the NBDC RDF portal side (as of 2017.11) so that it becomes an instance of a class representing the database to which the URI belongs.

The URI prefix of the main database to which RDF is published and the resource used when linking from the outside

Database name | Class | URI Prefix
------------- | ----- | ----------
uniprot | core: Protein | [http://purl.uniprot.org/uniprot/](http://purl.uniprot.org/uniprot/)
Ensembl | obo: SO_0001217 (protein_coding_gene) | [http://rdf.ebi.ac.uk/resource/ensembl/](http://rdf.ebi.ac.uk/resource/ensembl/)
ChEMBL | cco: Substance | [http://rdf.ebi.ac.uk/resource/chembl/molecule/](http://rdf.ebi.ac.uk/resource/chembl/molecule/)
ExpressionAtlas | atlas: BaseLineExpressionValue | [http://rdf.ebi.ac.uk/resource/expressionatlas/](http://rdf.ebi.ac.uk/resource/expressionatlas/)
 | atlas: DifferentialExpressionRatio | [http://rdf.ebi.ac.uk/resource/expressionatlas/](http://rdf.ebi.ac.uk/resource/expressionatlas/)
Reactome | biopax 3: Pathway | [http://identifiers.org/reactome/](http://identifiers.org/reactome/)
BioModels | | 	
BioSamples | biosd-terms: Sample | [http://rdf.ebi.ac.uk/resource/biosamples/sample](http://rdf.ebi.ac.uk/resource/biosamples/sample)
PubChem | compound | [http://rdf.ncbi.nlm.nih.gov/pubchem/compound/](http://rdf.ncbi.nlm.nih.gov/pubchem/compound/)
 | substance | [http://rdf.ncbi.nlm.nih.gov/pubchem/substance](http://rdf.ncbi.nlm.nih.gov/pubchem/substance)
MESH | meshv: Topical Descriptor | [http://id.nlm.nih.gov/mesh/](http://id.nlm.nih.gov/mesh/)
wwPDB | PDBo: datablock | [http://rdf.wwpdb.org/pdb/1NH2](http://id.nlm.nih.gov/mesh/)


##### 1.2.5 Add a link to document information

If you know the bibliographic information that serves as the basis for the content you describe, please actively attach a link to literature information. Link literature information using [PubMed](http://pubmed.org/) or [DOI ID](http://doi.org/) wherever possible . We recommend using the following prefix when URIizing these IDs.

	* URI of PubMed ID [http://rdf.ncbi.nlm.nih.gov/pubmed/24495517](http://rdf.ncbi.nlm.nih.gov/pubmed/24495517)
	* URI of DOI ID [http://doi.org/10.1021/jo 0349227](http://doi.org/10.1021/jo 0349227)


If the PubMed ID or if DOI ID is not available, Please describe as an instance of '''bibo:Article'''  of the academic literature [Bibliographic Ontology](http://bibliontology.com/). If it is a book, you can also use it '''bibo:Book''' . We also recommend using [Bibliographic Ontology](http://bibliontology.com/) [ [GitHub](https://github.com/structureddynamics/Bibliographic-Ontology-BIBO) ] for detailed information on literature (journal name, number of volumes, number of pages, publisher etc.).

When linking to literature information, we recommend using property [dcterms:references](http://purl.org/dc/terms/references) .

	'''<a resource> dcterms:references pubmed:24495517 .'''

We recommend using [Bibliographic Ontology](http://bibliontology.com/) for detailed information on the literature .

	'''
	@prefix dcterms: <http://purl.org/dc/terms/> .  
	@prefix bibo: <http://purl.org/ontology/bibo/> .  
	@prefix prism: <http://prismstandard.org/namespaces/1.2/basic/> .  

	<a resource> dcterms:references [  
  	 a bibo:Article;  
  	 prism:publicationName "Nature science cell";  
  	 prism:volume "10";  
  	 prism:number "11";  
  	 prism:startingPage "123";  
  	 prism:endingPage "456";  
  	 dcterms:date "2015-12-08" ;  
  	 seeAlso <http://rdf.ncbi.nlm.nih.gov/pubmed/1234567>  
	] .  '''

##### 1.2.6 Add meta information to data


One of the advantages of semanticizing the data is that you can give as much metadata as you need. By attaching metadata to the RDF data set, you can specify the creation date, person who created it, data source, category of data, licenses etc. There is [VoID](http://www.w3.org/TR/void/) as a vocabulary for describing metadata about RDF datasets . For example, [VoID Editor](http://voideditor.cs.man.ac.uk/) is a tool to efficiently construct metadata for RDF datasets using this vocabulary .

Also, for each statement (statement: each triple) described in RDF, you can give the provenance of that sentence. In particular, in the case of scientific data, since each sentence often expresses a fact, the history of writing that sentence is important. There are various possibilities in this history, from information described in papers and books, information obtained by mechanical method (homology search by BLAST etc.) to individual guessing. By explicitly describing the history, users will be able to select how to use each triple (such as using only the information mentioned in the paper).

Several vocabularies for describing log information in RDF have been proposed ( [DC terms](http://purl.org/dc/terms/) , [PROV-O](http://www.w3.org/TR/prov-o/) , [PAV](http://bioportal.bioontology.org/ontologies/PAV) ), and several models for describing log information have also been proposed ( [RDF Reification](https://goo.gl/RXvpBE) , [NanoPub](http://nanopub.org/wordpress/) , [OvoPub](http://arxiv.org/pdf/1305.6800.pdf) , [VoAG](http://linkedmodel.org/doc/voag/1.0/), etc.). BioHackathon 2014 was discussed as well.　 [Standardization of RDF data and development of tools / ontologies (page 8)](http://goo.gl/u32WCo).

##### 1.2.7 Link to images


* When [building](http://xmlns.com/foaf/spec/#term_depiction) a link to an image representing an entity that the subject URI means, we use [foaf:depiction](http://xmlns.com/foaf/spec/#term_depiction).

	'''an-assay-db:12345 foaf:depiction an-assay-db-image:12345.jpg .'''

* Conversely, if the URI of the image file is the subject and the resource URI of what is drawn in the image is the object, use [foaf:depicts](http://xmlns.com/foaf/spec/#term_depicts) .

	'''an-assay-db-image:12345.jpg foaf:depicts an-assay-db:12345 .'''

##### 1.2.8 Appropriate use of URI, blank node, literal properly

RDF consists of a combination of URI, blank node, and literal, but how to apply RDF conversion data to these will be examined to make it more useful RDF. The minimum unit of RDF is the subject, predicate, object triplet (RDF triple)

* Subject is URI or blank node
* Predicates include URI
* The object may be a URI, a blank node, or a literal

Each triple indicates that the subject and the object are connected by the relationship indicated by the predicate. Here, You can use,

* URI is appropriate when it is appropriate to globally uniquely identify the resources (mono or koto) indicated by the URI in the web world
* If it is not necessary to globally identify blank nodes such as merely combining specific RDF triples and joining them together. As an example, since UniProt requires that proteins identified by Q6GZX3 be globally identified, a URI of [http://purl.uniprot.org/uniprot/Q6GZX3](http://purl.uniprot.org/uniprot/Q6GZX3) is assigned.

Meanwhile, numerical data such as character strings and observation values ​​express the values ​​themselves, not identifiers (IDs), so they are expressed in literals. In this case, you can more clearly describe the meaning of the value (semantics of data) by attaching the data type such as unit to the value. In a string literal, you can specify languages ​​such as "English" and "Japanese" such as '''"protein"@en''' or '''"protein"@ja''' by using language tags, and to indicate that 123 is a numeric literal teeth '''"123"^^xsd:integer''' you can put the data type URI as (here '''xsd:integeris''' [http://www.w3.org/2001/XMLSchema#integer](http://www.w3.org/2001/XMLSchema#integer) of [QName](https://en.wikipedia.org/wiki/QName) is the notation). Furthermore, please refer to the "How to describe values ​​with units" section of this guideline for how to describe literals specifying units of numeric values.

##### 1.2.9 Do not use URIs starting with https.

Currently, many websites are progressing to HTTPS. This is a necessary effort to realize a secure WWW, but there is a problem with using a URI that starts with https when describing RDF. For example, if 　you use [http://identifiers.org/uniprot/Q6GZX3](http://identifiers.org/uniprot/Q6GZX3) for certain RDF data and [https://identifiers.org/uniprot/Q6GZX3](https://identifiers.org/uniprot/Q6GZX3)　for other RDF data as UniProt's Q6GZX3 URI , Although they describe the same resource, they are different URIs on the RDF data and can not be connected as it is. Building the RDF data while grasping all whether the specific website is https or whether the https or http scheme is used for each URI in the RDF data to be used, SPARQL Doing a search is not realistic. Therefore, we recommend that URIs beginning with 'http: //' be used as before for URIs used in RDF data.

####  1.3 Guidelines for using / building ontologies
##### 1.3.1 Reusing existing ontologies

When describing information as RDF, it is often not trivial that what ontology should be used for the subject class, predicate properties, etc. is a difficult problem. At least it is recommended to use it if there is something appropriate for the widely used ontology or vocabulary as follows.

**List of common vocabulary**
<table>
<thead>
<tr>
<th>vocabulary</th>
<th>Namespace</th>
<th>Reference link</th>
</tr>
</thead>
<tbody>
<tr>
<td>RDF 1.0</td>
<td><a href="http://www.w3.org/1999/02/22-rdf-syntax-ns#" rel="nofollow">http://www.w3.org/1999/02/22-rdf-syntax-ns#</a></td>
<td><a href="http://www.w3.org/TR/2004/REC-rdf-concepts-20040210/" rel="nofollow">concept and syntax</a> , [<a href="http://www.w3.org/TR/2004/REC-rdf-primer-20040210/" rel="nofollow">RDF Introduction </a>,  <a href="http://www.w3.org/TR/2004/REC-rdf-syntax-grammar-20040210/" rel="nofollow">RDF/XML</a></td>
</tr>
<tr>
<td>RDF 1.1</td>
<td><a href="http://www.w3.org/1999/02/22-rdf-syntax-ns#" rel="nofollow">http://www.w3.org/1999/02/22-rdf-syntax-ns#</a></td>
<td><a href="http://www.w3.org/TR/2014/REC-rdf11-concepts-20140225/" rel="nofollow">concept and syntax</a>, <a href="http://www.w3.org/TR/2014/NOTE-rdf11-primer-20140225/" rel="nofollow">RDF</a>, <a href="http://www.w3.org/TR/2014/REC-rdf-syntax-grammar-20140225/" rel="nofollow">RDF/XML</a>, <a href="http://www.w3.org/TR/2014/REC-turtle-20140225/" rel="nofollow">Turtle</a></td>
</tr>
<tr>
<td>RDFS 1.0</td>
<td><a href="http://www.w3.org/2000/01/rdf-schema#" rel="nofollow">http://www.w3.org/2000/01/rdf-schema#</a></td>
<td><a href="http://www.w3.org/TR/2004/REC-rdf-schema-20040210/" rel="nofollow">Specification</a> </td>
</tr>
<tr>
<td>RDFS 1.1</td>
<td><a href="http://www.w3.org/2000/01/rdf-schema#" rel="nofollow">http://www.w3.org/2000/01/rdf-schema#</a></td>
<td><a href="http://www.w3.org/TR/2014/REC-rdf-schema-20140225/" rel="nofollow">Specification</a></td>
</tr>
<tr>
<td>OWL 1</td>
<td><a href="http://www.w3.org/2002/07/owl" rel="nofollow">http://www.w3.org/2002/07/owl</a></td>
<td><a href="http://www.w3.org/TR/2004/REC-owl-features-20040210/" rel="nofollow">Summary</a> </td>
</tr>
<tr>
<td>OWL 2</td>
<td><a href="http://www.w3.org/2002/07/owl" rel="nofollow">http://www.w3.org/2002/07/owl</a></td>
<td><a href="http://www.w3.org/TR/2012/REC-owl2-overview-20121211/" rel="nofollow">Summary</a>]</td>
</tr>
<tr>
<td>DC</td>
<td><a href="http://purl.org/dc/elements/1.1/" rel="nofollow">http://purl.org/dc/elements/1.1/</a></td>
<td><a href="http://dublincore.org/documents/dces/" rel="nofollow">Specification</a></td>
</tr>
<tr>
<td>DC terms</td>
<td><a href="http://purl.org/dc/terms/" rel="nofollow">http://purl.org/dc/terms/</a></td>
<td> <a href="http://dublincore.org/documents/dcmi-terms/" rel="nofollow">Specification</a></td>
</tr>
<tr>
<td>SKOS</td>
<td><a href="http://www.w3.org/2004/02/skos/core" rel="nofollow">http://www.w3.org/2004/02/skos/core</a></td>
<td><a href="http://www.w3.org/TR/skos-reference/" rel="nofollow">Summary</a> , <a href="http://www.w3.org/TR/2009/NOTE-skos-primer-20090818/" rel="nofollow">Introduction to SKOS</a> </td>
</tr>
<tr>
<td>FOAF</td>
<td><a href="http://xmlns.com/foaf/0.1/" rel="nofollow">http://xmlns.com/foaf/0.1/</a></td>
<td><a href="http://xmlns.com/foaf/spec/" rel="nofollow">Specification</a>]</td>
</tr>
<tr>
<td>VoID</td>
<td><a href="http://rdfs.org/ns/void#" rel="nofollow">http://rdfs.org/ns/void#</a></td>
<td><a href="http://www.w3.org/TR/void/" rel="nofollow">Specification</a>,<a href="http://semanticweb.org/wiki/VoID" rel="nofollow">Summary</a>]</td>
</tr>
<tr>
<td>UO</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://code.google.com/p/unit-ontology/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/UO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>QUDT 1.1</td>
<td><a href="http://qudt.org/1.1/vocab/unit/" rel="nofollow">http://qudt.org/1.1/vocab/unit/</a></td>
<td><a href="http://www.linkedmodel.org/catalog/qudt/1.1/index.html" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/QUDT" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>QUDT 2.0</td>
<td><a href="http://qudt.org/2.0/schema/qudt/" rel="nofollow">http://qudt.org/2.0/schema/qudt/</a></td>
<td><a href="http://qudt.org/doc/2017/DOC_SCHEMA-QUDT-v2.0.html" rel="nofollow">Home</a></td>
</tr>
<tr>
<td>PROV-O</td>
<td><a href="http://www.w3.org/ns/prov#" rel="nofollow">http://www.w3.org/ns/prov#</a></td>
<td><a href="http://www.w3.org/TR/prov-o/" rel="nofollow">Specification</a>, <a href="http://bioportal.bioontology.org/ontologies/PROVO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>PAV</td>
<td><a href="http://purl.org/pav/" rel="nofollow">http://purl.org/pav/</a></td>
<td><a href="http://www.essepuntato.it/lode/http://purl.org/pav/2.0/" rel="nofollow">Specification</a>], <a href="http://bioportal.bioontology.org/ontologies/PAV" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>XSD</td>
<td><a href="http://www.w3.org/2001/XMLSchema#" rel="nofollow">http://www.w3.org/2001/XMLSchema#</a></td>
<td></td>
</tr>
<tr>
<td>DCAT</td>
<td><a href="http://www.w3.org/ns/dcat#" rel="nofollow">http://www.w3.org/ns/dcat#</a></td>
<td><a href="https://www.w3.org/TR/vocab-dcat/" rel="nofollow">Specification</a></td>
</tr>
<tr>
<td>BIBO</td>
<td><a href="http://purl.org/ontology/bibo/" rel="nofollow">http://purl.org/ontology/bibo/</a></td>
<td><a href="http://bibliontology.com/" rel="nofollow">Home</a>, <a href="https://github.com/structureddynamics/Bibliographic-Ontology-BIBO">GitHub</a></td>
</tr>
<tr>
<td>Event</td>
<td><a href="http://purl.org/NET/c4dm/event.owl#" rel="nofollow">http://purl.org/NET/c4dm/event.owl#</a></td>
<td><a href="http://motools.sourceforge.net/event/event.html" rel="nofollow">Home</a></td>
</tr>
<tr>
<td>GEO</td>
<td><a href="http://www.w3.org/2003/01/geo/wgs84_pos#" rel="nofollow">http://www.w3.org/2003/01/geo/wgs84_pos#</a></td>
<td><a href="https://www.w3.org/2003/01/geo/" rel="nofollow">Home</a></td>
</tr>
</tbody>
</table>

[http://www.linkedmodel.org/catalog/qudt/1.1/index.html](http://www.linkedmodel.org/catalog/qudt/1.1/index.html)
In addition, you can find commonly used vocabulary in [Linked Open Vocabularies (LOV)](http://lov.okfn.org/dataset/lov/) .

In addition to ontologies and vocabularies that handle general concepts as described above, there are also ontologies that are easy to use when describing life science information. Especially [BioPortal](http://bioportal.bioontology.org/) can search suitable ontology and vocabulary in life science.

**Domain ontology on life science information**
<table>
<thead>
<tr>
<th>abbreviation</th>
<th>Ontology name</th>
<th>Namespace</th>
<th>Reference link</th>
</tr>
</thead>
<tbody>
<tr>
<td>GO</td>
<td>Gene Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://geneontology.org/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/GO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>PRO</td>
<td>Protein Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://pir.georgetown.edu/pro/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/PR" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>SO</td>
<td>Sequence Types and Features Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://www.sequenceontology.org/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/SO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>FALDO</td>
<td>Feature Annotation Location Description Ontology</td>
<td><a href="http://biohackathon.org/resource/faldo#" rel="nofollow">http://biohackathon.org/resource/faldo#</a></td>
<td><a href="https://github.com/JervenBolleman/FALDO">GitHub</a></td>
</tr>
<tr>
<td>PO</td>
<td>Plant Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://plantontology.org/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/PO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>Taxonomy</td>
<td>Taxnomy Ontology</td>
<td><a href="http://ddbj.nig.ac.jp/ontologies/taxonomy/" rel="nofollow">http://ddbj.nig.ac.jp/ontologies/taxonomy/</a></td>
<td><a href="http://bioportal.bioontology.org/ontologies/NCBITAXON" rel="nofollow">BioPortal</a>, <a href="http://tga.nig.ac.jp/ontologies/" rel="nofollow">DDBJ</a></td>
</tr>
<tr>
<td>Nucleotide</td>
<td>INSDC Nucleotide Sequence Entry Ontology</td>
<td><a href="http://ddbj.nig.ac.jp/ontologies/nucleotide/" rel="nofollow">http://ddbj.nig.ac.jp/ontologies/nucleotide/</a></td>
<td><a href="http://tga.nig.ac.jp/ontologies/" rel="nofollow">DDBJ</a></td>
</tr>
<tr>
<td>MeSH</td>
<td>Medical Subject Headings</td>
<td></td>
<td><a href="http://bioportal.bioontology.org/ontologies/MESH" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>CL</td>
<td>Cell Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://bioportal.bioontology.org/ontologies/CL" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>FMA</td>
<td>Foundational Model of Anatomy</td>
<td><a href="http://purl.org/sig/ont/fma/" rel="nofollow">http://purl.org/sig/ont/fma/</a></td>
<td><a href="http://sig.biostr.washington.edu/projects/fm/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/FMA" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>UBERON</td>
<td>Uber Anatomy Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://uberon.org" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/UBERON" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>SNOMED-CT</td>
<td>Systematized Nomenclature of Medicine - Clinical Terms</td>
<td><a href="http://purl.bioontology.org/ontology/SNOMEDCT/" rel="nofollow">http://purl.bioontology.org/ontology/SNOMEDCT/</a></td>
<td><a href="http://bioportal.bioontology.org/ontologies/SNOMEDCT" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>SIO</td>
<td>Semanticscience Integrated Ontology</td>
<td><a href="http://semanticscience.org/resource/" rel="nofollow">http://semanticscience.org/resource/</a></td>
<td><a href="http://semanticscience.org/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/SIO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>EFO</td>
<td>Experimental Factor Ontology</td>
<td><a href="http://www.ebi.ac.uk/efo/" rel="nofollow">http://www.ebi.ac.uk/efo/</a></td>
<td><a href="http://www.ebi.ac.uk/efo" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/EFO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>ECO</td>
<td>Evidence Ontology</td>
<td><a href="http://purl.obolibrary.org/obo" rel="nofollow">http://purl.obolibrary.org/obo</a></td>
<td><a href="http://code.google.com/p/evidenceontology/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/ECO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>EDAM</td>
<td>EDAM bioinformatics operations, data types, formats, identifiers and topics</td>
<td><a href="http://edamontology.org/" rel="nofollow">http://edamontology.org/</a></td>
<td><a href="http://edamontology.org/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/EDAM" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>CMO</td>
<td>Clinical Measurement Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://phenoonto.sourceforge.net/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/CMO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>MMO</td>
<td>Measurement Method Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://phenoonto.sourceforge.net/" rel="nofollow">Home</a>, <a href="http://purl.bioontology.org/ontology/MMO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>XCO</td>
<td>Experimental Conditions Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://phenoonto.sourceforge.net/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/XCO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>OrthO</td>
<td>Ortholog Ontology</td>
<td><a href="http://purl.jp/bio/11/orth#" rel="nofollow">http://purl.jp/bio/11/orth#</a></td>
<td><a href="http://mbgd.genome.ad.jp/ontology/" rel="nofollow">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/ORTHO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>PIERO</td>
<td>PIERO Enzyme Reaction Ontology　</td>
<td></td>
<td><a href="http://reactionontology.org/" rel="nofollow">Home</a></td>
</tr>
<tr>
<td>GlycoRDF</td>
<td>Glycan Ontology</td>
<td><a href="http://purl.jp/bio/12/glyco/glycan#" rel="nofollow">http://purl.jp/bio/12/glyco/glycan#</a></td>
<td><a href="https://github.com/ReneRanzinger/GlycoRDF">Home</a>, <a href="http://bioportal.bioontology.org/ontologies/GLYCORDF" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>MONDO</td>
<td>Monarch Disease Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="http://www.obofoundry.org/ontology/mondo.html" rel="nofollow">Home</a>, <a href="https://bioportal.bioontology.org/ontologies/MONDO" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>HPO</td>
<td>The Human Phenotype Ontology</td>
<td><a href="http://purl.obolibrary.org/obo/" rel="nofollow">http://purl.obolibrary.org/obo/</a></td>
<td><a href="https://hpo.jax.org/" rel="nofollow">Home</a>, <a href="https://bioportal.bioontology.org/ontologies/HP" rel="nofollow">BioPortal</a></td>
</tr>
<tr>
<td>HCO</td>
<td>The Human Chromosome Ontology</td>
<td><a href="http://identifiers.org/hco/" rel="nofollow">http://identifiers.org/hco/</a></td>
<td><a href="https://github.com/med2rdf/hco">github</a></td>
</tr>
</tbody>
</table>

As a vocabulary concerning the language name, [ISO 639-1](http://id.loc.gov/vocabulary/iso639-1.html) or [ISO 639-2](http://id.loc.gov/vocabulary/iso639-2.html) provided by the National Diet Library can be used.

When using URIs of existing ontologies and vocabularies in SPARQL and Turtle, we often write abbreviated namespaces. For example, the above '''rdf:typeURI''' is '''<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>''' a shorthand notation. This abbreviated form is written according to the notation method of Turtle's prefixed name, and if correspondence between abbreviated name and actual notation is contradictory declared in one dataset, any short name can be used Although there is no http://prefix.cc/abbreviated name widely used for major ontologies, using it makes the human readable. You can use [prefix.cc]() as a service to look for general abbreviated names .


##### 1.3.2 Building a new ontology

Ontology is sometimes required when converting data and knowledge into RDF. Typical examples are:

* When explicitly indicating the class to which the resource belongs
* Although some information is described as a literal, if it is likely to be conceptualized
* When a predicate that describes the relationship between subject resource and object resource is needed

Etc. From the viewpoint of data integration, it is preferable to use an existing ontology, but if you can not find an appropriate one among the existing ones, you only have to build a new ontology. On the semantic web, we generally define an ontology using an ontology description language called OWL. We also write OWL itself in RDF. You can write OWL (RDF too) with a text editor, but you can work more efficiently by using a dedicated ontology editor. In particular, an ontology editor which can do work with GUI is indispensable when constructing an ontology . The following editors can describe the OWL ontology.

	* [Protege](http://protege.stanford.edu/) Open source ontology editor
	* [WebProtege](http://webprotege.stanford.edu/)	cloud version of Protege
	* [TopBraid](http://www.topquadrant.com/tools/modeling-topbraid-composer-standard-edition/) Commercial ontology editor

##### 1.3.3 Defining domain and range for properties

For properties (predicate of RDF; predicate), '''rdfs:domain''' and '''rdfs:range''' can be defined. This makes the meaning of the property more explicit.

	'''example）  
	core:classifiedWith rdfs:domain core:Protein .  
	core:classifiedWIth rdfs:range core:Concept . ''' 

In the above example, '''core:classifiedWith''' property gives functional annotation information such as GO terms to the protein entry of UniProt . it has a domain '''core:Protein''' and '''core:Concept''' is defined as a range . From the string "classifiedWith", humans can pick up certain meanings such as "classify subjects by something", but machines do not do that.

In the machine,

	'''uniprot:P51028 core:classifiedWith go:2000044 . '''

Also,
	
	'''uniprot:P51028 ex:fugahoge go:2000044 .  '''

Also, that between '''uniprot:P51028''' and the '''go:2000044''' the statement does not have the meaning of the above relationship, .

By domain and range defined in the ontology, this property '''core:classifiedWIth''' is  a '''core:Protein''' not taken subject to resources other than type, '''core:Concept''' that can not have a value other than the machine is able to be used in the information processing.

##### 1.3.4 Properly describe ontology classes and properties

When new ontology classes are defined, it is important to give appropriate explanation so that the user can understand what the concept represented by the class is. As a property, please write clearly because it is simple and easy to use '''skos:definition''', '''rdfs:comment''' etc.

	'''example)  
	core:Active_Site_Annotation rdfs:comment "Amino acid(s) involved in the activity of an enzyme."^^xsd:string; '''
	
#### 1.4 Guidelines for Providing RDF
##### 1.4.1 Regular release and version information

Currently, I think that many databases are converting from databases that have already been constructed and provided in another format (relational database etc.) when converting to RDF. In such a case, there will be a problem of synchronization on the contents between the data provided as a public service and the data converted to RDF. As much as possible, periodic RDF conversion is desirable so that the difference between the data provided to the primary and the RDF converted data does not become large. For that reason, please convert from existing format to RDF conversion as much as possible, and from that process eliminate as much as possible the steps that require manual work. When a large amount of manual work occurs in the process of RDF conversion, periodic RDF conversion becomes extremely difficult. Manual work (annotation, ontology mapping etc) should be included in the primary database construction work. Of course, this is not the case if the primary database format is RDF.

Also, when converting to RDF, it is necessary to give an appropriate version number. For versioning , you can use pav:version (for RDF data) or owl:versionInfo(for ontology).

	* [pav: version](http://purl.org/pav/version)
	* [owl: versionInfo](http://www.w3.org/2002/07/owl#versionInfo)


##### 1.4.2 Add license information

By explicitly writing the license, the user can use the RDF data with confidence. From the viewpoint of ease of use of data, it is recommended that you use an appropriate [license](https://creativecommons.org/licenses/) of Creative Commons . Although there is room for discussion about ontology originally whether it is a work or not, assuming the case of copyrighted work, and from the ease of reusing a part, if CC 0 is good It is said.


##### 1.4.3 To provide a schema diagram
##### 1.4.4 Provide SPARQL sample
##### 1.4.5 CORS countermeasure when publishing SPARQL endpoint

#### 1.5 RDF model easy to use in life science
##### 1.5.1 Method of describing measurement items and measurement values
##### 1.5.2 How to describe gene and protein sequence coordinate information
##### 1.5.3 How to write samples
