"""
Project 2 — Biological Data Analysis

Implement all functions below. Do not change their signatures.
"""

import json
import csv
import math
import re
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Task 1 — JSON: Genomic Variants (ClinVar)
# ─────────────────────────────────────────────────────────────────────────────

CHROM_LENGTHS = {
    "1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,
    "6":170805979,"7":159345973,"8":145138636,"9":138394717,"10":133797422,
    "11":135086622,"12":133275309,"13":114364328,"14":107043718,"15":101991189,
    "16":90338345,"17":83257441,"18":80373285,"19":58617616,"20":64444167,
    "21":46709983,"22":50818468
}

# nao podemos usar counter??? PODEMOS ! TO DO MUDAR!

def mostPathogenicGene(filepath):
    symbols= []
    tmp={}
    with open(filepath, 'r') as f:
        data= list(json.load(f)["result"].values())
    data= data[1:]
    for i in range(len(data)):
        entry = dict(data[i])
        # it returns gene as type list but knows how to read it as a dict
        # $$ genes é uma lista com um dictionario dentro $$
        if entry["clinical_significance"]["description"] == "Pathogenic":
            # nao me parece que counter seja aplicavel atras disto abaixo (devido a $$)
            for gene in entry["genes"]: #dict itrs gene
                tmp[gene["symbol"]] = tmp.get(gene["symbol"], 0) + 1 # nao tiver default é por 0
    counter =Counter(tmp)
    return (counter.most_common()[0])
    #---------OU----------
    #gene_symbol, count = "", 0
    #if tmp[entry["genes"]["symbol"]] >= count:
    #   gene_symbol, count = entry["genes"]["symbol"], tmp[entry["genes"]["symbol"]]
    #return gene_symbol, count
#print(mostPathogenicGene("data\clinvar_variants.json"))


# verifiquei que sao todos: "assembly": "GRCh38"
# a identificação do autosome é chr, igues...
'''
most are like this:
    "chr": "17",
              "start": 43166073,
              "stop": 43166073,
              "assembly": "GRCh38
so the variant only happens in one specific base...?
              '''
def variantDensityByChrom(filepath):
    # GRCh38 chromosome lengths (base pairs)
    CHROM_LENGTHS = {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "8": 145138636,
        "9": 138394717,
        "10": 133797422,
        "11": 135086622,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "19": 58617616,
        "20": 64444167,
        "21": 46709983,
        "22": 50818468
    }
    counts = {}
    with open(filepath, 'r') as f:
        data= json.load(f)["result"]
    # de certeza que ha outra maneira sem ser bues for loops
    for uid in data["uids"]:
        entry = data[uid]
        for variation in entry["variation_set"]:
            for loc in variation["variation_loc"]:
                chrom = loc["chr"]
                if  loc["assembly"] == "GRCh38" and chrom in CHROM_LENGTHS:
                        counts[chrom] = counts.get(chrom, 0) + 1 # count variants by chromossomes
                        break # se tiver a mesma variaçao, mais do que uma vez registada
    densities = {
    chrom: round(count / (CHROM_LENGTHS[chrom] / 1_000_000), 4)
    for chrom, count in counts.items()}
    return dict(sorted(densities.items(), key=lambda x: x[1], reverse=True))
    # lambda because the default its sorting by key

#density = variantDensityByChrom('data/clinvar_variants.json')
#for chrom, d in list(density.items())[:5]:
#    print(f"  chr{chrom}: {d} variants/Mbp")

def topPathogenicConditions(filepath):
    """T1.3 — Returns [(condition, count)] for conditions with >=6 Pathogenic/Likely pathogenic variants."""
    with open(filepath, 'r') as f:
        data= json.load(f)["result"]
    counts ={}
    for uid in data["uids"]:
        entry = data[uid]
        significance = entry["clinical_significance"]["description"] 
        # assign significance istead of for loop, to be able to if-condictioning the significance and still be at entry["condictions"] level!
        if significance in ("Pathogenic", "Likely pathogenic"):
            for condition in entry["conditions"]:
                counts[condition["name"]] = counts.get(condition["name"], 0) +1 
        #make sure >= 6:
    res = [(name, count) for name, count in counts.items() if count >= 6]
        #descending order by value count (same last ex)
    return sorted(res, key=lambda x: x[1], reverse= True)

for cond, n in topPathogenicConditions('data/clinvar_variants.json')[:5]:
    print(f"  {cond}: {n}")

def variantIndex(filepath):
    """T1.4 — Returns nested dict {gene: {chr: {clinical_significance: [titles]}}}."""
    with open(filepath, 'r') as f:
        data= json.load(f)["result"]
    res = {}
    for uid in data["uids"]:
        entry = data[uid]
        title = entry["title"]
        significance = entry["clinical_significance"]["description"]
        for gene in entry["genes"]:
            gene_symbol= gene["symbol"]
            for variation in entry["variation_set"]:
                for loc in variation["variation_loc"]:
                    chromosome = loc["chr"]
                    #acho que é mesmo preciso gerar cada dictionario de cada vez
                    if gene_symbol not in res:
                        res[gene_symbol] = {}
                    if chromosome not in res[gene_symbol]:
                        res[gene_symbol][chromosome] = {}
                    if significance not in res[gene_symbol][chromosome]:
                        res[gene_symbol][chromosome][significance] = []
                    res[gene_symbol][chromosome][significance].append(title)
    return res
    
# level1 : extract gene; variation_set; title
# aka; gene["symbol"]:{ at same level of gene) variation_set[{variation_loc}] <- usar for: loc["chr"] : { clinical_significance["name"] : title}}}
# That is: `{ gene_symbol : { chromosome : { clinical_significance : [variant_titles] } } }`
idx = variantIndex('data/clinvar_variants.json')
for gene in list(idx.keys())[:2]:
    for chrom in idx[gene]:
        for sig, titles in idx[gene][chrom].items():
            print(f"  {gene} | chr{chrom} | {sig}: {len(titles)} variant(s)")
print("------------------")
print(list(idx.items())[:5]        )
print(len(idx['BRCA2']['13']['Pathogenic']) )

##TO DO generalize in aux functions (like open_json + skip top levels that are irrelevant; get_gene: get_chr...) 

# ─────────────────────────────────────────────────────────────────────────────
# Task 2 — NumPy: Gene Expression Analysis (TCGA)
# ─────────────────────────────────────────────────────────────────────────────

def geneExpressionStats(expr_matrix):
    """T2.1 — Returns (means, stds): two 1D arrays of shape (n_genes,)."""
    pass

def topDifferentialGenes(expr_matrix, gene_ids, tumor_mask, normal_mask):
    """T2.2 — Returns list of 10 tuples (gene_id, delta_mean) sorted by descending |delta_mean|."""
    pass

def sampleCorrelationMatrix(expr_matrix):
    """T2.3 — Returns Pearson correlation matrix of shape (n_samples, n_samples). NumPy only."""
    pass

def zscoreNormalize(expr_matrix):
    """T2.4a — Returns z-score normalised matrix (per gene, across all samples)."""
    pass

def classifyGenes(z_matrix, gene_ids, tumor_mask):
    """T2.4b — Returns {gene_id: 'overexpressed'|'underexpressed'|'stable'}."""
    pass
