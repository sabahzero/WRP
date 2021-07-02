# python3 WRP.py to run in terminal
## runs

# Nodes
import pandas

from pathlib import Path
from tqdm.autonotebook import tqdm 

from data_tools.df_processing import char_combine_iter 
from data_tools.wiki import node_query_pipeline


nodes = []


# Anatomy 
q = """SELECT DISTINCT ?anatomy ?anatomyLabel 
        WHERE {
          ?anatomy wdt:P1554 ?uberon
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }""" 

res = node_query_pipeline(q, {}, 'anatomy')
nodes.append(res)

# Biological Process 
q = """SELECT DISTINCT ?biological_process ?biological_processLabel 
        WHERE {
          ?biological_process wdt:P31 wd:Q2996394 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'biological_process')
nodes.append(res)

# Cellular Component
q = """SELECT DISTINCT ?cellular_component ?cellular_componentLabel 
    WHERE {
      ?cellular_component wdt:P31 wd:Q5058355 .
      SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
    }"""

res = node_query_pipeline(q, {}, 'cellular_component')
nodes.append(res)

# Compounds 
q = """SELECT DISTINCT ?compound ?compoundLabel
        WHERE {
          ?compound wdt:P31 wd:Q11173 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }
        limit 150000""" 

res = node_query_pipeline(q, {}, 'compound')
nodes.append(res)

# Disease
q = """SELECT DISTINCT ?disease ?diseaseLabel 
        WHERE {
          ?disease wdt:P31 wd:Q12136 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""
 
res = node_query_pipeline(q, {}, 'disease')
nodes.append(res)

# Genes (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?gene ?geneLabel
        WHERE {
          ?gene wdt:P31 wd:Q7187 .
          ?gene wdt:P703 wd:Q15978631 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'gene')
nodes.append(res)

# Pathway
q = """SELECT DISTINCT ?pathway ?pathwayLabel
        WHERE {
          ?pathway wdt:P31 wd:Q4915012 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'pathway')
nodes.append(res)

# Phenotype (nothing for hpo? apply to Compound?)
q = """SELECT DISTINCT ?phenotype ?phenotypeLabel ?hpo 
        WHERE {
          {?phenotype wdt:P31 wd:Q169872.}UNION{?phenotype wdt:P3841 ?hpo}
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'phenotype')
nodes.append(res)

# Protein (note focus on Homo sapiens) 
q = """SELECT DISTINCT ?protein ?proteinLabel
        WHERE {
          ?protein wdt:P31 wd:Q8054 .
          ?protein wdt:P703 wd:Q15978631 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'protein')
nodes.append(res)

# Molecular Function
q = """SELECT DISTINCT ?molecular_function ?molecular_functionLabel 
        WHERE {
          ?molecular_function wdt:P31 wd:Q14860489 .
          SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }
        }"""

res = node_query_pipeline(q, {}, 'molecular_function')
nodes.append(res)


nodes = pandas.concat(nodes, sort=False, ignore_index=True)

out_dir = Path('../results/')
out_dir.mkdir(parents=True, exist_ok=True)
nodes.to_csv(out_dir.joinpath('01a_nodes.csv'), index=False)