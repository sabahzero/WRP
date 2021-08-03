# Run in terminal with command 'python3 WRP.py' *in src folder
## Make sure you've first run "pip install -r 00_requirements.txt"

import pandas as pd

import time # sleep function
import functools # what does this do?
from pathlib import Path
from itertools import chain # what does this do?
from tqdm.autonotebook import tqdm 

from data_tools.df_processing import char_combine_iter, add_curi
from data_tools.wiki import execute_sparql_query, node_query_pipeline, standardize_nodes, standardize_edges


# Nodes (works - Jul 13, 2021)
nodes = []

## Anatomy 
q = """SELECT DISTINCT ?anatomy ?anatomyLabel 
        WHERE {
          ?anatomy wdt:P1554 ?uberon
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }""" 

res = node_query_pipeline(q, {}, 'anatomy')
nodes.append(res)

## Biological Process 
q = """SELECT DISTINCT ?biological_process ?biological_processLabel 
        WHERE {
          ?biological_process wdt:P31 wd:Q2996394 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'biological_process')
nodes.append(res)

## Cellular Component
q = """SELECT DISTINCT ?cellular_component ?cellular_componentLabel 
    WHERE {
      ?cellular_component wdt:P31 wd:Q5058355 .
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
    }"""

res = node_query_pipeline(q, {}, 'cellular_component')
nodes.append(res)

## Compounds 
q = """SELECT DISTINCT ?compound ?compoundLabel
        WHERE {
          ?compound wdt:P31 wd:Q11173 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }
        limit 150000""" 

res = node_query_pipeline(q, {}, 'compound')
nodes.append(res)

## Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel 
        WHERE {
          ?disease wdt:P31 wd:Q12136 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""
 
res = node_query_pipeline(q, {}, 'disease')
nodes.append(res)

## Genes (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?gene ?geneLabel
        WHERE {
          ?gene wdt:P31 wd:Q7187 .
          ?gene wdt:P703 wd:Q15978631 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'gene')
nodes.append(res)

## Pathway
q = """SELECT DISTINCT ?pathway ?pathwayLabel
        WHERE {
          ?pathway wdt:P31 wd:Q4915012 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'pathway')
nodes.append(res)

## Phenotype (note focus on Homo sapiens)
q = """SELECT DISTINCT ?phenotype ?phenotypeLabel 
        WHERE {
          {?phenotype wdt:P31 wd:Q169872.}UNION{?phenotype wdt:P3841 ?hpo}
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'phenotype')
nodes.append(res)

## Protein (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?protein ?proteinLabel
        WHERE {
          ?protein wdt:P31 wd:Q8054 .
          ?protein wdt:P703 wd:Q15978631 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'protein')
nodes.append(res)

## Molecular Function
q = """SELECT DISTINCT ?molecular_function ?molecular_functionLabel 
        WHERE {
          ?molecular_function wdt:P31 wd:Q14860489 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'molecular_function')
nodes.append(res)

## Compile all 10 node type categories
nodes = pd.concat(nodes, sort=False, ignore_index=True)

out_dir = Path('../results/')
out_dir.mkdir(parents=True, exist_ok=True)
nodes.to_csv(out_dir.joinpath('01a_nodes.csv'), index=False)

# Edges
## Put in a sleep function (for 5 seconds) to override the SPARQL overflow

def process_taxa(edges): # Integrate process_taxa() function into data_tools package ?
    nodes = edges.drop_duplicates(subset=['taxon', 'tax_id'])[['taxon', 'taxonLabel', 'tax_id']]
    nodes = add_curi(nodes, {'tax_id': 'NCBITaxon'})
    return standardize_nodes(nodes, 'taxon')

# What is happening in this code cell?
# Why do we need nodes to get edges? Is it a good idea that we have them?

prev_dir = Path('../results/').resolve()
prev_nodes = pd.read_csv(prev_dir.joinpath('01a_nodes.csv')) 

nodes = []
edges = []

# Approach 1
## Syntax 1 -- Direct statement: Disease causes infection
q = """SELECT DISTINCT ?disease ?taxon ?taxonLabel ?tax_id
    WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
      ?disease p:P828 [ps:P828 wd:Q166231;pq:P642 ?taxon;].
      OPTIONAL{?taxon wdt:P685 ?tax_id}.
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}"""

qr = execute_sparql_query(q) # Query
tax_nodes = process_taxa(qr) # Query by taxa
edge_res = standardize_edges(qr, 'taxon', 'disease', 'causes') # Standardize taxon, disease, and causes
nodes.append(tax_nodes) # Update nodes
edges.append(edge_res) # Update edges

## Syntax 2 -- Qualifier statements
### a. disease has-cause TAXON 
q = """SELECT DISTINCT ?disease ?diseaseLabel ?doid ?taxon ?taxonLabel ?tax_id
    WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
        ?taxon wdt:P685 ?tax_id. 
       {?disease wdt:P828 ?taxon}UNION{?taxon wdt:P1542 ?disease}.
        OPTIONAL {?disease wdt:P699 ?doid.}
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}"""

qr = execute_sparql_query(q)
tax_nodes = process_taxa(qr)
edge_res = standardize_edges(qr, 'taxon', 'disease', 'causes')
nodes.append(tax_nodes)
edges.append(edge_res)

### b. TAXON has-effect Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel ?doid ?taxon ?taxonLabel ?tax_id
    WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
        ?taxon wdt:P685 ?tax_id.
           {?disease wdt:P828 ?taxon}UNION{?taxon wdt:P1542 ?disease}.
           OPTIONAL {?disease wdt:P699 ?doid.}
           SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}"""

qr = execute_sparql_query(q)
tax_nodes = process_taxa(qr)
edge_res = standardize_edges(qr, 'taxon', 'disease', 'causes')
nodes.append(tax_nodes)
edges.append(edge_res)

# Approach 2
## Syntax 1
q = """SELECT DISTINCT ?disease ?diseaseLabel ?doid ?parent_tax ?parent_taxLabel ?par_taxid ?taxon ?taxonLabel ?tax_id
    WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
      ?disease p:P828 [ps:P828 wd:Q166231;
                       pq:P642 ?parent_tax;].
      OPTIONAL{?disease wdt:P699 ?doid}.
      OPTIONAL{?parent_tax wdt:P685 ?par_taxid}.
      FILTER NOT EXISTS {?parent_tax wdt:P105 wd:Q36732}.
      FILTER NOT EXISTS {?parent_tax wdt:P105 wd:Q3978005}.
      {?taxon wdt:P171+ ?parent_tax}UNION{?parent_tax wdt:P171+ ?taxon}
      ?taxon wdt:P105 wd:Q7432 .
      ?taxon wdt:P685 ?tax_id    
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}"""

qr = execute_sparql_query(q)
tax_nodes = process_taxa(qr)
edge_res = standardize_edges(qr, 'taxon', 'disease', 'causes', 'computed')
edge_res['comp_type'] = 'punning' # What does this do? What is meant by 'punning'?
nodes.append(tax_nodes)
edges.append(edge_res)

## Syntax 2 
q = """SELECT DISTINCT ?disease ?diseaseLabel ?doid ?parent_tax ?parent_taxLabel ?parent_tax_id ?taxon ?taxonLabel ?tax_id
    WHERE 
    {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
        ?parent_tax wdt:P685 ?parent_tax_id. 
      FILTER NOT EXISTS {?parent_tax wdt:P105 wd:Q36732}.
      FILTER NOT EXISTS {?parent_tax wdt:P105 wd:Q3978005}.      
       {?disease wdt:P828 ?parent_tax}UNION{?parent_tax wdt:P1542 ?disease}.
        OPTIONAL {?disease wdt:P699 ?doid.}
      {?taxon wdt:P171+ ?parent_tax}UNION{?parent_tax wdt:P171+ ?taxon}
      ?taxon wdt:P685 ?tax_id .
      ?taxon wdt:P105 wd:Q7432 .
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}"""

qr = execute_sparql_query(q)
tax_nodes = process_taxa(qr)
edge_res = standardize_edges(qr, 'taxon', 'disease', 'causes', 'computed')
edge_res['comp_type'] = 'punning' 
nodes.append(tax_nodes)
edges.append(edge_res)

# Remove duplicates
tax_nodes = pd.concat(nodes, sort=False, ignore_index=True).drop_duplicates(subset=['id'])
nodes = [tax_nodes]

# IN TAXON (needs to come first for ENCODES)
## Focuses on taxa with annotations to genes or proteins in Wikidata
### Genes 
q = """SELECT DISTINCT ?taxon
    WHERE {?gene wdt:P31 wd:Q7187.
      ?gene wdt:P703 ?taxon.}"""

qr = execute_sparql_query(q)
gene_taxa = set(qr['taxon'])

q = """SELECT DISTINCT ?gene ?geneLabel ?entrez ?symbol ?hgnc ?omim ?ensembl
        WHERE {{
          ?gene wdt:P31 wd:Q7187.
          ?gene wdt:P703 wd:{tax}.
          OPTIONAL{{?gene wdt:P351 ?entrez .}}
          OPTIONAL{{?gene wdt:P353 ?symbol .}}
          OPTIONAL{{?gene wdt:P354 ?hgnc .}}
          OPTIONAL{{?gene wdt:P492 ?omim .}}
          OPTIONAL{{?gene wdt:P594 ?ensembl .}}
          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}}}"""

tax_gene_edges = []
gene_curi_map = {'entrez': 'NCBIGene', 'symbol': 'SYM', 'hgnc':'HGNC', 'omim':'OMIM', 'ensembl':'ENSG'}

for tax_id in gene_taxa & set(tax_nodes['id']):
    this_q = q.format(tax=tax_id)
    res = node_query_pipeline(this_q, gene_curi_map, 'gene')
    if res is None:
        continue
    nodes.append(res[['id', 'name', 'label', 'xrefs']].copy())
    res['tax'] = tax_id
    res_edges = standardize_edges(res, 'id', 'tax', 'in_taxon')
    tax_gene_edges.append(res_edges)
    
gene_tax = pd.concat(tax_gene_edges, sort=False, ignore_index=True)
edges.append(gene_tax)

q = """SELECT DISTINCT ?taxon
    WHERE {?gene wdt:P31 wd:Q7187.
      ?gene wdt:P703 ?taxon.}"""

qr = execute_sparql_query(q)
gene_taxa = set(qr['taxon'])

q = """SELECT DISTINCT ?gene ?geneLabel ?entrez ?symbol ?hgnc ?omim ?ensembl
        WHERE {{
          ?gene wdt:P31 wd:Q7187.
          ?gene wdt:P703 wd:{tax}.
          OPTIONAL{{?gene wdt:P351 ?entrez .}}
          OPTIONAL{{?gene wdt:P353 ?symbol .}}
          OPTIONAL{{?gene wdt:P354 ?hgnc .}}
          OPTIONAL{{?gene wdt:P492 ?omim .}}
          OPTIONAL{{?gene wdt:P594 ?ensembl .}}
          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}}}"""

tax_gene_edges = []
gene_curi_map = {'entrez': 'NCBIGene', 'symbol': 'SYM', 'hgnc':'HGNC', 'omim':'OMIM', 'ensembl':'ENSG'}

for tax_id in gene_taxa & set(tax_nodes['id']):
    this_q = q.format(tax=tax_id)
    res = node_query_pipeline(this_q, gene_curi_map, 'gene')
    if res is None:
        continue
    nodes.append(res[['id', 'name', 'label', 'xrefs']].copy())
    res['tax'] = tax_id
    res_edges = standardize_edges(res, 'id', 'tax', 'in_taxon')
    tax_gene_edges.append(res_edges)
    
gene_tax = pd.concat(tax_gene_edges, sort=False, ignore_index=True)
edges.append(gene_tax)

### Proteins
q = """SELECT DISTINCT ?taxon
    WHERE {?protein wdt:P31 wd:Q8054.
      ?protein wdt:P703 ?taxon.}"""

qr = execute_sparql_query(q)
prot_taxa = set(qr['taxon'])

q = """SELECT DISTINCT ?protein ?proteinLabel ?uniprot
        WHERE {{
          ?protein wdt:P31 wd:Q8054.
          ?protein wdt:P703 wd:{tax}.
          OPTIONAL{{?protein wdt:P352 ?uniprot .}}
          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}}}"""

tax_prot_edges = []

for tax_id in prot_taxa & set(tax_nodes['id']):
    this_q = q.format(tax=tax_id)
    res = node_query_pipeline(this_q, {'uniprot':'UniProt'}, 'protein')
    if res is None:
        continue
    nodes.append(res[['id', 'name', 'label', 'xrefs']].copy())
    res['tax'] = tax_id
    res_edges = standardize_edges(res, 'id', 'tax', 'in_taxon')
    tax_prot_edges.append(res_edges)
    
prot_tax = pd.concat(tax_prot_edges, sort=False, ignore_index=True)
edges.append(prot_tax)

# ASSOCIATED WITH
## Gene ASSOCIATED WITH Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel ?gene ?geneLabel 
WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
    ?gene wdt:P31 wd:Q7187 .
    ?disease wdt:P2293 ?gene
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'gene', 'disease', 'associated_with')
edges.append(edge_res)

## Pathway ASSOCIATED WITH Disease
q = """SELECT DISTINCT ?pathway ?disease 
WHERE {
    ?pathway wdt:P31 wd:Q4915012 .
    FILTER NOT EXISTS{?pathway wdt:P686 ?goid}
    {?disease wdt:P31 wd:Q12136 .}UNION{?disease wdt:P699 ?doid}
    {?pathway wdt:P1050 ?disease}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'pathway', 'disease', 'associated_with')
edges.append(edge_res)

# ENABLES
## Protein ENABLES Molecular Function
q = """SELECT DISTINCT ?protein ?molecular_function
WHERE {?protein wdt:P31 wd:Q8054.
  ?molecular_function wdt:P686 ?goid.
  {?protein wdt:P680 ?molecular_function}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'protein', 'molecular_function', 'enables')
edges.append(edge_res)

# ENCODES
## Gene ENCODES Protein (not focus on Homo sapiens)
q = """SELECT DISTINCT ?gene ?protein 
WHERE {{
  ?gene wdt:P31 wd:Q7187.
  ?gene wdt:P703 wd:{tax}.
  ?protein wdt:P31 wd:Q8054.
  ?protein wdt:P703 wd:{tax}.
  {{?gene wdt:P688 ?protein}}UNION{{?protein wdt:P702 ?gene}}}}"""

human_tax_id = 'Q15978631'
encodes_edges = []
infectious_tax = list(set(gene_tax['end_id']) & set(prot_tax['end_id']))

for tax in infectious_tax + [human_tax_id]:
    this_q = q.format(tax=tax)
    qr = execute_sparql_query(this_q)
    if qr is not None:
        this_edge = standardize_edges(qr, 'gene', 'protein', 'encodes')
        encodes_edges.append(this_edge)
        
encodes_edges = pd.concat(encodes_edges, sort=False, ignore_index=True)
edges.append(encodes_edges) # Why encodes_edges vs edge_res?

# HAS PART
## Pathway HAS PART Compoound
q = """SELECT DISTINCT ?pathway ?compound 
WHERE {
    ?pathway wdt:P31 wd:Q4915012 .
    FILTER NOT EXISTS{?pathway wdt:P686 ?goid}
    ?compound wdt:P31 wd:Q11173 .
    ?pathway wdt:P527 ?compound}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'pathway', 'compound', 'has_part')
edges.append(edge_res)

## Pathway HAS PART Gene
q = """SELECT DISTINCT ?pathway ?gene 
WHERE {?pathway wdt:P31 wd:Q4915012 .
    FILTER NOT EXISTS{?pathway wdt:P686 ?goid}
    ?gene wdt:P31 wd:Q7187 .
    ?pathway wdt:P527 ?gene}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'pathway', 'gene', 'has_part')
edges.append(edge_res)

# INVOLVED IN
## Pathway INVOLVED IN Biological Process
q = """SELECT DISTINCT ?pathway ?bio_process 
WHERE {
    ?pathway wdt:P31 wd:Q4915012 .
    FILTER NOT EXISTS{?pathway wdt:P686 ?goid}
    ?bio_process wdt:P31 wd:Q2996394 .
    ?bio_process wdt:P686 ?goid .
    {?pathway wdt:P31 ?bio_process}UNION{?bio_process wdt:P31 ?pathway}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'pathway', 'bio_process', 'involved_in')
edges.append(edge_res)

## Protein INVOLVED IN Biological Process
q = """SELECT DISTINCT ?protein ?biological_process
WHERE {?protein wdt:P31 wd:Q8054.
  ?biological_process wdt:P686 ?goid.
  {?protein wdt:P682 ?biological_process}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'protein', 'biological_process', 'involved_in')
edges.append(edge_res)

# INTERACTS WITH
## Compound INTERACTS WITH Protein
q = """SELECT DISTINCT ?compound ?compoundLabel ?qualifier ?qualifierLabel ?protein ?proteinLabel 
WHERE {
    ?compound wdt:P31 wd:Q11173 .
    ?protein wdt:P31 wd:Q8054 .
    { ?compound p:P129 [ps:P129 ?protein;
                  pq:P2868 ?qualifier] }
    UNION { ?protein p:P129 [ps:P129 ?compound;
                 pq:P366 ?qualifier] }
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'compound', 'protein', 'qualifierLabel')
edges.append(edge_res)

# PART OF
## Protein PART OF Cellular Component
q = """SELECT DISTINCT ?protein ?cell_component
WHERE {?protein wdt:P31 wd:Q8054.
  ?cell_component wdt:P686 ?goid.
  {?protein wdt:P681 ?cell_component}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'protein', 'cell_component', 'part_of')
edges.append(edge_res)

# SITE OF
## Anatomy SITE OF Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel ?anatomy ?anatomyLabel 
WHERE {
    {?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
    {?disease1 wdt:P31 wd:Q12136}UNION{?disease1 wdt:P699 ?doid1}.
    ?anatomy wdt:P1554 ?uberon .
    {?disease1 wdt:P927 ?anatomy} {?disease wdt:P279? ?disease1}
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'anatomy', 'disease', 'site_of')
edges.append(edge_res)

# ? Presents ? is this an edge..?
## Disease to Phenotype

q = """SELECT DISTINCT ?disease ?pheno
  WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}
      {?pheno wdt:P31 wd:Q169872.}UNION{?pheno wdt:P3841 ?hpo}
      {?pheno wdt:P780 ?disease}UNION{?disease wdt:P780 ?pheno}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'disease', 'pheno', 'presents')
edges.append(edge_res)

# ?? SUBCLASS OF ?? What is the purpose of this 'later punning' meaning...?
## Anatomy SUBCLASS OF Anatomy
q = """SELECT DISTINCT ?anatomy ?anatomyLabel ?anatomy1 ?anatomy1Label
WHERE {
    ?anatomy wdt:P1554 ?uberon .
    ?anatomy1 wdt:P1554 ?uberon1 .
    ?anatomy wdt:P279? ?anatomy1
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'anatomy', 'anatomy1', 'subclass_of')
edges.append(edge_res)

## Disease SUBCLASS OF Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel ?disease1 ?disease1Label
WHERE {{?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
    {?disease1 wdt:P31 wd:Q12136}UNION{?disease1 wdt:P699 ?doid1}.
    ?disease wdt:P279? ?disease1
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'disease', 'disease1', 'subclass_of')
edges.append(edge_res)

# TREATS
## Compound TREATS Disease
q = """SELECT DISTINCT ?compound ?compoundLabel ?disease ?diseaseLabel
  WHERE {
    ?compound wdt:P31 wd:Q11173 .
    {?disease wdt:P31 wd:Q12136}UNION{?disease wdt:P699 ?doid}.
    {?compound wdt:P2175 ?disease}UNION{?disease wdt:P2176 ?compound}
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'compound', 'disease', 'treats')
edges.append(edge_res)

## Compound TREATS Phenotype
q = """SELECT DISTINCT ?compound ?pheno
  WHERE {?compound wdt:P31 wd:Q11173 .
      {?pheno wdt:P31 wd:Q169872.}UNION{?pheno wdt:P3841 ?hpo}
      {?pheno wdt:P2176 ?compound}UNION{?compound wdt:P2175 ?pheno}}"""

qr = execute_sparql_query(q)
edge_res = standardize_edges(qr, 'compound', 'pheno', 'treats')
edges.append(edge_res)

nodes.append(prev_nodes)

nodes = pd.concat(nodes, sort=False, ignore_index=True)
edges = pd.concat(edges, sort=False, ignore_index=True)

out_dir = Path('../results/')
out_dir.mkdir(parents=True, exist_ok=True)

nodes.to_csv(out_dir.joinpath('01b_nodes.csv'), index=False)
edges.to_csv(out_dir.joinpath('01b_edges.csv'), index=False)
