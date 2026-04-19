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

#because i repeat this step em cada uma das tasks:
def load_variants(filepath):
    with open(filepath, 'r') as f:
        result = json.load(f)["result"]
    return result["uids"], result

CHROM_LENGTHS = {
    "1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,
    "6":170805979,"7":159345973,"8":145138636,"9":138394717,"10":133797422,
    "11":135086622,"12":133275309,"13":114364328,"14":107043718,"15":101991189,
    "16":90338345,"17":83257441,"18":80373285,"19":58617616,"20":64444167,
    "21":46709983,"22":50818468
}


def mostPathogenicGene(filepath):

    uids, result = load_variants(filepath)
    counter = Counter()
    for uid in uids:
        entry = result[uid]
        # genes é uma lista com um dictionario dentro
        if entry["clinical_significance"]["description"] == "Pathogenic":
            # nao me parece que counter seja aplicavel atras disto abaixo (devido a $$)
            for gene in entry["genes"]: #dict itrs gene
                counter[gene["symbol"]] +=1
    return (counter.most_common(1)[0])

    #---------OU----------
    #gene_symbol, count = "", 0
    #if tmp[entry["genes"]["symbol"]] >= count:
    #   gene_symbol, count = entry["genes"]["symbol"], tmp[entry["genes"]["symbol"]]
    #return gene_symbol, count

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

def variantDensityByChrom(filepath):
    # GRCh38 chromosome lengths (base pairs)
    
    counts = Counter()
    uids, result = load_variants(filepath)
    for uid in uids:
        entry = result[uid]
        for variation in entry["variation_set"]:
            for loc in variation["variation_loc"]:
                chrom = loc["chr"]
                if  loc["assembly"] == "GRCh38" and chrom in CHROM_LENGTHS:
                        counts[chrom] += 1
                        break # se tiver a mesma variaçao, nao ser mais do que uma vez registada
    densities = {
    chrom: round(count / (CHROM_LENGTHS[chrom] / 1_000_000), 4)
    for chrom, count in counts.items()}

    return dict(sorted(densities.items(), key=lambda x: x[1], reverse=True))
    # lambda because the default its sorting by key


def topPathogenicConditions(filepath):
    """T1.3 — Returns [(condition, count)] for conditions with >=6 Pathogenic/Likely pathogenic variants."""
    uids, result= load_variants(filepath)
    counts = Counter()
    for uid in uids:
        entry = result[uid]
        significance = entry["clinical_significance"]["description"] #nao era necessario quebrar mas para ser mais facil visualizar
        if significance in ("Pathogenic", "Likely pathogenic"):
            for condition in entry["conditions"]:
                counts[condition["name"]] += 1 #isto nao necessario chamar logo counter primeiro : counts.get(condition["name"], 0) +1 
        #make sure >= 6:
    res = [(name, count) for name, count in counts.items() if count >= 6]
        #descending order by value count (same last ex)
    return sorted(res, key=lambda x: x[1], reverse= True)


def variantIndex(filepath):
    """T1.4 — Returns nested dict {gene: {chr: {clinical_significance: [titles]}}}."""
    uids, result = load_variants(filepath)
    idx = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for uid in uids:
        entry = result[uid]
        title = entry["title"]
        significance = entry["clinical_significance"]["description"]
        chroms = set()
        for variation in entry["variation_set"]:
            for loc in variation["variation_loc"]:
                if loc["assembly"] == "GRCh38":
                    chroms.add(loc["chr"])
        for gene in entry["genes"]:
            gene_symbol = gene["symbol"]
            for chrom in chroms:
                idx[gene_symbol][chrom][significance].append(title)
    res ={}
    for gene, chroms in idx.items():
        res[gene] = {}
        for chrom, sig in chroms.items():
            res[gene][chrom] = dict(sig)
    return res
    
# level1 : extract gene; variation_set; title
# aka; gene["symbol"]:{ at same level of gene) variation_set[{variation_loc}] <- usar for: loc["chr"] : { clinical_significance["name"] : title}}}
# That is: `{ gene_symbol : { chromosome : { clinical_significance : [variant_titles] } } }`


## DONE TO DO generalize in aux functions (like open_json + "skip" top levels that are irrelevant; get_gene: get_chr...) 

# ─────────────────────────────────────────────────────────────────────────────
# Task 2 — NumPy: Gene Expression Analysis (TCGA)
# ─────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt
# Read expression matrix


with open('data/expression_matrix.csv', newline='') as f:
    reader = csv.reader(f)
    header = next(reader)
    sample_ids = header[1:]            # list of sample IDs
    gene_ids   = []
    expr_rows  = []
    for row in reader:
        gene_ids.append(row[0])
        expr_rows.append([float(x) for x in row[1:]])

expr_matrix = np.array(expr_rows)     # shape: (n_genes, n_samples)
gene_ids    = np.array(gene_ids)
sample_ids  = np.array(sample_ids)

# Read sample metadata
import pandas as pd
metadata = pd.read_csv('data/sample_metadata.csv')

# Boolean masks for tumor and normal samples
tumor_mask  = np.array([s in metadata[metadata['condition']=='tumor']['sample_id'].values
                         for s in sample_ids])
normal_mask = np.array([s in metadata[metadata['condition']=='normal']['sample_id'].values
                         for s in sample_ids])

def geneExpressionStats(expr_matrix):
    """T2.1 — Returns (means, stds): two 1D arrays of shape (n_genes,)."""
    # the standard qoulde give the results overall(all genes and all smaples)
    # by columns (axis=0) this is samples
    means = np.mean(expr_matrix, axis=1)
    stds = np.std(expr_matrix, axis=1)
    return means, stds


def topDifferentialGenes(expr_matrix, gene_ids, tumor_mask, normal_mask):
    """T2.2 — Returns list of 10 tuples (gene_id, delta_mean) sorted by descending |delta_mean|."""
    # mean expression by tumour SAMPLES (colums)
    mean_tumor = np.mean(expr_matrix[:, tumor_mask], axis=1)
    mean_normal = np.mean(expr_matrix[:, normal_mask], axis=1)
    delta_mean = mean_tumor - mean_normal

    delta_idx= np.argsort(abs(delta_mean))[::-1][:10]

    '''[:10] substitui e depois return adpatado
    lst= []
    for i in delta_idx[:10]:
        lst.append((gene_ids[i], delta_mean[i]))'''
    return [(gene_ids[i], delta_mean[i]) for i in delta_idx]


def sampleCorrelationMatrix(expr_matrix):
    """T2.3 — Returns Pearson correlation matrix of shape (n_samples, n_samples). NumPy only."""
    #matrix normal:(n_genes, n_samples); transposta: (n_samples, n_genes)
    #give the (n_samples, n_samples), correlation of rows!
    return np.corrcoef(expr_matrix.T)


def zscoreNormalize(expr_matrix):
    """T2.4a — Returns z-score normalised matrix (per gene, across all samples)."""
    means = np.mean(expr_matrix, axis=1, keepdims= True)
    stds= np.std(expr_matrix, axis=1, keepdims= True)

    z = np.zeros_like(expr_matrix, dtype=float)
    no_zero = (stds[:, 0] != 0)
    z[no_zero] = (expr_matrix[no_zero] -means[no_zero]) / stds[no_zero]
    return z



def classifyGenes(z_matrix, gene_ids, tumor_mask):
    """T2.4b — Returns {gene_id: 'overexpressed'|'underexpressed'|'stable'}."""
    genes= {}
    tumor_means = np.mean(z_matrix[:, tumor_mask], axis = 1)
    for i, gene in enumerate(gene_ids):
        if tumor_means[i] > 0.5:
            genes[gene] = 'overexpressed'
        elif tumor_means[i] < -0.5:
            genes[gene] = 'underexpressed'
        else:
            genes[gene] = 'stable'
    return genes



#TESTS

#---T.1
#1.1
print("mostPathogenicGene", mostPathogenicGene('data/clinvar_variants.json'))
print("mostPathogenicGene type", type(mostPathogenicGene('data\clinvar_variants.json')))

#1.2
density = variantDensityByChrom('data/clinvar_variants.json')
print("Density: top 5")
for chrom, d in list(density.items())[:5]:
    print(f"  chr{chrom}: {d} variants/Mbp")
print("type return", type(density))

#1.3
print("top Pathogenic condictions! top 5")
for cond, n in topPathogenicConditions('data/clinvar_variants.json')[:5]:
    print(f"  {cond}: {n}")

print("rtuenr type", type(topPathogenicConditions('data/clinvar_variants.json')))

#1.4
print("variant Index")
idx = variantIndex('data/clinvar_variants.json')
for gene in list(idx.keys())[:2]:
    for chrom in idx[gene]:
        for sig, titles in idx[gene][chrom].items():
            print(f"  {gene} | chr{chrom} | {sig}: {len(titles)} variant(s)")
print("------------------")
print(list(idx.items())[:5]        )
print(len(idx['BRCA2']['13']['Pathogenic']) )

#---T.2

print(f"Expression matrix shape: {expr_matrix.shape}")
print(f"Tumor samples: {tumor_mask.sum()} | Normal samples: {normal_mask.sum()}")

#2.1
print("gene expr stats")
means, stds = geneExpressionStats(expr_matrix)
print(means.shape, stds.shape)

means, stds = geneExpressionStats(expr_matrix)
print("Top 5 most variable genes (highest std):")
top_var_idx = np.argsort(stds)[::-1][:5]
for i in top_var_idx:
    print(f" GENE_EXPRESSION\n {gene_ids[i]}: mean={means[i]:.3f}, std={stds[i]:.3f}")

#2.2
top10 = topDifferentialGenes(expr_matrix, gene_ids, tumor_mask, normal_mask)
print("Top 10 differentially expressed genes:")
for gene, delta in top10:
    direction = "▲ overexpressed" if delta > 0 else "▼ underexpressed"
    print(f"  {gene}: Δmean = {delta:+.3f}  ({direction} in tumor)")

print(top10[0]) # not the same of the example but em baxio vemos que esta top

#2.3
corr = sampleCorrelationMatrix(expr_matrix)
print(f"Correlation matrix shape: {corr.shape}")
print(f"Min correlation: {corr.min():.4f} | Max (off-diagonal): {corr[corr < 0.9999].max():.4f}")
# Visualise as heatmap — complete the plot below

fig, ax = plt.subplots(figsize=(9, 7))
sh= ax.imshow(corr, cmap='coolwarm', vmin=-1, vmax=1)
n_tumor = tumor_mask.sum()
ax.axvline(n_tumor - 0.5, color='black', lw=1)
ax.axhline(n_tumor - 0.5, color='black', lw=1)
ax.set_title("Pearson Correlation (tumour vs normal)")

plt.colorbar(sh, ax=ax)
plt.tight_layout()
plt.savefig('sample_correlation_heatmap.png', dpi=150)
plt.show()

#2.4
#zscore
z= zscoreNormalize(expr_matrix)
print("Zscore:\n", z)
#last
classes = classifyGenes(z, gene_ids, tumor_mask)

counts = Counter(classes.values())
print("Gene classification summary:")
for label, n in sorted(counts.items()):
    print(f"  {label}: {n} gene(s)")
print("\nDifferentially expressed genes (overexpressed or underexpressed):")
for gene, cls in classes.items():
    #if cls != 'stable':
        print(f"  {gene}: {cls}")