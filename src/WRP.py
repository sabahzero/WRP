# Run in terminal with command 'python3 WRP.py' *in src folder *use venv as needed
## Make sure you've first run "pip install -r 00_requirements.txt"

# Relevant Packages

import pandas as pd

import time 
from datetime import datetime

from pathlib import Path
from tqdm.autonotebook import tqdm 

from data_tools.df_processing import char_combine_iter 
from data_tools.wiki import node_query_pipeline

# Create time stamp of run
timeStringNow = datetime.now().strftime("+%Y-%m-%dT00:00:00Z") 
start_time = time.time()

# Nodes (works - Aug 5, 2021)
nodes = []

## Anatomy 
q = """SELECT DISTINCT ?anatomy ?anatomyLabel ?uberon ?mesh
WHERE {
  ?anatomy wdt:P1554 ?uberon  
           OPTIONAL{?anatomy wdt:P486 ?mesh .}
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}""" 

res = node_query_pipeline(q, {'uberon':'UBERON', 'mesh': 'MESH'}, 'anatomy')
nodes.append(res)

## Biological Process 
q = """SELECT DISTINCT ?biological_process ?biological_processLabel ?goid
WHERE {
  ?biological_process wdt:P31 wd:Q2996394 .
  ?biological_process wdt:P686 ?goid
                      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'goid':'GO'}, 'biological_process')
nodes.append(res)

## Cellular Component
q = """SELECT DISTINCT ?cellular_component ?cellular_componentLabel ?goid
WHERE {
  ?cellular_component wdt:P31 wd:Q5058355 .
  ?cellular_component wdt:P686 ?goid
                      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'goid':'GO'}, 'cellular_component')
nodes.append(res)

## Compounds 
q = """SELECT DISTINCT ?compound ?compoundLabel ?kegg_drug ?chebi ?drugbank_id ?umlscui ?chembl_id ?unii ?ikey ?pubchem_cid ?rxnorm ?mesh_supplemental_record_ui ?mesh_descriptor_ui
WHERE {
  ?compound wdt:P31 wd:Q11173 .
  OPTIONAL { ?compound wdt:P665 ?kegg_drug .}
  OPTIONAL { ?compound wdt:P683 ?chebi .}
  OPTIONAL { ?compound wdt:P715 ?drugbank_id .}
  OPTIONAL { ?compound wdt:P2892 ?umlscui .}
  OPTIONAL { ?compound wdt:P592 ?chembl_id .}
  OPTIONAL { ?compound wdt:P652 ?unii .}
  OPTIONAL { ?compound wdt:P3350 ?ikey .}
  OPTIONAL { ?compound wdt:P662 ?pubchem_cid .}
  OPTIONAL { ?compound wdt:P3345 ?rxnorm .}
  OPTIONAL { ?compound wdt:P6680 ?mesh_supplemental_record_ui .}
  OPTIONAL { ?compound wdt:P486 ?mesh_descriptor_ui .}
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}
limit 17000""" 

res = node_query_pipeline(q, {'unii': 'UNII', 'rxnorm': 'RxCUI', 'drugbank_id': 'DB', 
                              'umlscui': 'UMLS', 'chebi': 'CHEBI', 'chembl_id': 'CHEMBL',
                              'kegg_drug': 'KEGG', 'ikey': 'IKEY', 'pubchem_cid': 'PCID', 
                              'mesh_supplemental_record_ui': 'MESH', 
                              'mesh_descriptor_ui': 'MESH'}, 'compound')
nodes.append(res)

## Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel ?umlscui ?snomed_ct ?doid ?mesh ?mondo ?omim ?orpha
WHERE {{
  ?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
       OPTIONAL {?disease wdt:P2892 ?umlscui .}
       OPTIONAL {?disease wdt:P5806 ?snomed_ct. }
       OPTIONAL {?disease wdt:P699 ?doid. }
       OPTIONAL {?disease wdt:P486 ?mesh. }
       OPTIONAL {?disease wdt:P5270 ?mondo. }
       OPTIONAL {?disease wdt:P492 ?omim. }
       OPTIONAL {?disease wdt:P1550 ?orpha. }
       SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""
    
res = node_query_pipeline(q, {'umlscui': 'UMLS', 'snomed_ct': 'SNOMED', 'mesh': 'MESH',
                              'doid': 'DOID', 'mondo': 'MONDO', 'omim': 'OMIM', 
                              'orpha': 'ORPHA'}, 'disease')
nodes.append(res)

## Genes (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?gene ?geneLabel ?entrez ?symbol ?hgnc ?omim ?ensembl
WHERE {
  ?gene wdt:P31 wd:Q7187 .
  ?gene wdt:P703 wd:Q15978631 .
  OPTIONAL{{?gene wdt:P351 ?entrez .}}
  OPTIONAL{{?gene wdt:P353 ?symbol .}}
  OPTIONAL{{?gene wdt:P354 ?hgnc .}}
  OPTIONAL{{?gene wdt:P492 ?omim .}}
  OPTIONAL{{?gene wdt:P594 ?ensembl .}}
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'entrez': 'NCBIGene', 'symbol': 'SYM', 'hgnc':'HGNC', 
                              'omim':'OMIM', 'ensembl':'ENSG'}, 'gene')
nodes.append(res)

## Pathway
q = """SELECT DISTINCT ?pathway ?pathwayLabel ?react ?wpid
WHERE {
  ?pathway wdt:P31 wd:Q4915012 .
  OPTIONAL{?pathway wdt:P3937 ?react .}
  OPTIONAL{?pathway wdt:P2410 ?wpid .}
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'react':'REACT', 'wpid':'WP'}, 'pathway')
nodes.append(res)

## Phenotype (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?phenotype ?phenotypeLabel ?hpo ?mesh ?omim ?snomed  
WHERE {{
  ?phenotype wdt:P31 wd:Q169872.}UNION{?phenotype wdt:P3841 ?hpo}
       OPTIONAL {?phenotype wdt:P3841 ?hpo .}
       OPTIONAL {?phenotype wdt:P486 ?mesh . }
       OPTIONAL {?phenotype wdt:P492 ?omim . }
       OPTIONAL {?phenotype wdt:P5806 ?snomed . }
       SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'mesh': 'MESH', 'omim': 'OMIM', 
                              'hpo':'HP', 'snomed': 'SNOMED'}, 'phenotype')
nodes.append(res)

## Protein (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?protein ?proteinLabel ?uniprot
WHERE {
  ?protein wdt:P31 wd:Q8054 .
  ?protein wdt:P703 wd:Q15978631 .
  OPTIONAL{{?protein wdt:P352 ?uniprot .}}
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'uniprot':'UniProt'}, 'protein')
nodes.append(res)

## Molecular Function
q = """SELECT DISTINCT ?molecular_function ?molecular_functionLabel ?goid
WHERE {
  ?molecular_function wdt:P31 wd:Q14860489 .
  ?molecular_function wdt:P686 ?goid
                      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

res = node_query_pipeline(q, {'goid':'GO'}, 'molecular_function')
nodes.append(res)

## Concatenate and convert to .csv
nodes = pd.concat(nodes, sort=False, ignore_index=True)

out_dir = Path('../results/')
out_dir.mkdir(parents=True, exist_ok=True)

nodes.to_csv(out_dir.joinpath('01a_nodes.csv'), index=False)

# Output and print when query is complete
end_time = time.time() 
print("The total time of this query is:", (end_time - start_time)/60, "minutes")