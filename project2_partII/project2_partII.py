"""
Project 2 — Biological Data Analysis (Part II)

Implement all functions below. Do not change their signatures.
"""

import json
import csv
import math
import re
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import networkx as nx


# ─────────────────────────────────────────────────────────────────────────────
# Task 3 — pandas - GWAS Catalogue
# ─────────────────────────────────────────────────────────────────────────────

def mostStudiedTraits(gwas):
    """
    T3.1 — Returns a DataFrame with the traits with >= 3 studies, sorted by descending number of studies.
    """
    # TODO
    pass

def mostSignificantPerChrom(gwas):
    """
    T3.2 — Returns a dictionary with the most significant SNP on each chromosome (ie., with lowest p-value).
    """
    # TODO
    pass

def publicationTrend(gwas):
    """
    T3.3 — Returns a DataFrame indexed by year with columns ['n_studies', 'n_unique_traits', 'n_unique_snps'], 
    restricted to years with >= 5 studies. Sorted by year.
    """
    # TODO
    pass

"""
T3.3 — Plot n_studies over time
"""
# TODO

def genomicHotspots(gwas):
    """
    T3.4 — Returns the top-5 1-Mbp windows with the most distinct associated traits, as a list of tuples 
    (chr, window_start, window_end, n_traits).
    """
    # TODO
    pass


# ─────────────────────────────────────────────────────────────────────────────
# Task 4 — NetworkX - Protein–Protein Interaction Network
# ─────────────────────────────────────────────────────────────────────────────

#tsv its like csv but istead of commas as tab separating
def df_tsv(filepath):
    df = pd.read_csv(filepath, sep='\t')
    # print(df.columns) <- its right
    return df
#print(df_tsv("data/ppi_network.tsv"))

def buildPPINetwork(filepath, min_score = 400):
    #check one by one combined_score >= min_score
    #OR redefine df
    df = df_tsv(filepath) # generalizaçao bue minima
    df = df[df['combined_score'] >= min_score] #needs the outer df[] inside its just :
    '''
0       True
1       True
2       True
3       True
4       True
    ...  
223    False
'''
#each node protein1 or protein2 
# edge combined_score
    Gr = nx.Graph()
    Gr.add_weighted_edges_from(df[['protein_a', 'protein_b', 'combined_score']].values)
    return Gr
G= buildPPINetwork("data/ppi_network.tsv")
print((G.number_of_nodes()), G.number_of_edges())

def topHubs(graph):
    """
    T4.1 — Returns a list of the 10 proteins with the highest degree, sorted by descending degree.
    """
    degrees = graph.degree
    sort = sorted(degrees, key=lambda x: x[1], reverse=True)
    return sort[:10]

hubs = topHubs(G)
print("Top 10 network hubs:")
for protein, deg in hubs:
    print(f"  {protein}: {deg} interaction partners")


def networkComponents(graph):
    """
    T4.2 — Returns (n_components, largest_size, component_sizes_sorted_desc).
    """
    # componenents aka vertices...
    # count vertices...
    components = sorted(nx.connected_components(graph), key= len, reverse= True)
    component_sizes_sorted_desc = [len(comp) for comp in components]
    largest_size = component_sizes_sorted_desc[0]
   # n_components = nx.number_connected_components(graph)
    n_components = len(components)
    #print(type(n_components)) #isto dava jeito de visualizar....
    #print(component_sizes_sorted_desc)
    return (n_components, largest_size, component_sizes_sorted_desc)

print(networkComponents(G))


def shortestInteractionPath(graph, protein_a, protein_b):
    """
    T4.3 — Returns the list of proteins along the shortest (unweighted) path between protein_a and protein_b. 
    Returns None if no path exists.
    """
    path =nx.shortest_path(graph, protein_a, protein_b)
    return path

path = (shortestInteractionPath(G, "TP53", "EGFR"))
if path:
    print(" → ".join(path))
    print(f"Path length: {len(path)-1} edge(s)")
else:
    print("No path found.")
print(shortestInteractionPath(G, "TP53", "APC"))
print(shortestInteractionPath(G, "BIRC2", "APC"))
print(shortestInteractionPath(G,  "APC", "BIRC2"))


def breastCancerModule(graph, clinvar_filepath):
    #ler json: gerar list, dict, set: only variants Pathogenic and likely pathogenic
    # SUB GRAGH
    
    pass

# Test
subG, mod = breastCancerModule(G, 'data/clinvar_variants.json')
print(f"Breast cancer module: {subG.number_of_nodes()} proteins, {subG.number_of_edges()} interactions")
if mod is not None:
    print(f"Modularity score: {mod:.4f}")
else:
    print("Module too small for community detection.")


# Visualise the subgraph
fig, ax = plt.subplots(figsize=(8, 6))
pos = nx.spring_layout(subG, seed=42)
nx.draw_networkx(subG, pos=pos, ax=ax,
                 node_color='salmon', node_size=800,
                 font_size=9, edge_color='gray', width=1.5)
ax.set_title("Breast Cancer Disease Module (PPI Network)", fontsize=12)
ax.axis('off')
plt.tight_layout()
plt.savefig('breast_cancer_module.png', dpi=150)
plt.show()