import re
from collections import Counter
#import Biopython #just to use for checks

## Part I

### T1

def readFASTA(file: str):
    with open(file, 'r',encoding='utf-8-sig') as f:
        d = {}
        ID = None
        seq = ""
        desc= ""
        for line in f:
            line = line.rstrip("\n")
            if line.startswith('>'):
                if ID is not None:
                    d[ID] = (desc, seq)
                header = line[1:].strip()
                spl= header.split(maxsplit=1)
                ID = spl[0]
                if len(spl) > 1: desc = spl[1]
                else: desc=""
                seq = ""
            else:
                seq += line
            if ID is not None:
                    d[ID] = (desc, seq)
    return d

### T2

def seqComposition(seq, type):
    if type == "DNA":
        dict = {'A':0, 'C':0, 'G':0, 'T':0} # ordem de exemplo
    
    elif type == "RNA":
        dict = {'A':0, 'C':0, 'G':0, 'U': 0}
    
    elif type == "PROTEIN":
        dict = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0,
                'E':0, 'Q': 0, 'H': 0, 'I':0, 'L': 0,
                'K':0, 'M': 0, 'F': 0, 'P': 0, 'S':0,
                'T':0, 'W': 0, 'Y':0, 'V': 0}
    else:
        return None
    for nuc in seq:
        if nuc in dict.keys():
            dict[nuc] += 1
    
    return dict

def dnaGCcontent(seq):
   dict = seqComposition(seq, "DNA")
   cg = dict['C'] + dict['G']
   total = dict['C'] + dict['G'] + dict['T'] + dict['A'] # tal como seqcomp ignorar chars maus
   return str(round((cg/total * 100), 2)) + ('%'+ " GC composition") 

## Parte II

### T3

def transcribeDNA2RNA(seq):
    rna_seq = ""
    for nuc in seq:
        if nuc == 'T':
            rna_seq+='U'
        else:
            rna_seq+=nuc
    return rna_seq     

### T4

def translateRNA2protein(seq):
    genetic_code =[]
    dict ={}
    protein =""
    with open("data_sequences/genetic_code.txt", 'r', encoding='utf-8-sig') as f:
        for line in  f:
            line.rstrip("\n")
            temp =line.split()
            temp2 = temp[-1]
            genetic_code.append(temp2)
        AAs = genetic_code[0]
        base_1 = genetic_code[2]
        base_2 = genetic_code[3]
        base_3 = genetic_code[4]
        for i in range(0, len(base_1)-1):
            dict[base_1[i]+base_2[i]+base_3[i]] = AAs[i]
        for j in range(0, len(seq)-2, 3):
            codon = str(seq[j] + seq[j+1] + seq[j+2])
            if codon in dict.keys():
                protein += dict[codon]
        return protein

### T5

def reverseComplement(sequence : str):
    dna_complement = {'A': 'T', 'C':'G', 'T': 'A', 'G':'C'}
    dna_com = []
    for base in sequence:
        if base in dna_complement.keys():
            dna_com.append(dna_complement[base])
    return "".join(dna_com[::-1])

## Parte III

### T6

def findMotif(dna_sequence : str, motif :str):
    step = len(motif)
    positions = []
    for i in range(len(dna_sequence)-step + 1):
        if dna_sequence[i:i+step] == motif:
            positions.append(i)

    return positions

### T7


def mostFrequentKMotifs(dna_sequence: str, k : int):
    motifs = []
    result = []
    for b in range(0, len(dna_sequence) -k +1):
            motifs.append(dna_sequence[b: b + k])
    if not motifs:
        return []
    most_motifs = Counter(motifs)
    max_freq = most_motifs.most_common(1)[0][1]
    for i,y in most_motifs.items():
        if y == max_freq:
            result.append((i,y))
    if len(result) == 1:
        return result[0]
    return result

## Parte IV

### T8

def highestGC(file):
    fasta = readFASTA(file)
    max_gc = -1
    max_id = None
    for id, (_, seq) in fasta.items():
        comp = seqComposition(seq, "DNA")
        cg = comp['C'] + comp['G']
        total = sum(comp.values())
        if total == 0:
            continue
        gc_percent = round((cg / total) * 100)
        if gc_percent > max_gc:
            max_gc = gc_percent
            max_id = id
    return (max_id, max_gc)

### T9

def compositionMatrix(file):
    fasta = readFASTA(file)
    result = {}
    for id, (_, seq) in fasta.items():
        comp = seqComposition(seq, "DNA")
        result[id] = comp
    return result

### T10

def indexSpeciesGC(fasta: dict):
    out_dict = {}

    for id, (description, seq) in fasta.items():
        found = re.search(r"\[(.*?)\]", description)
        if found:
            specie = found.group(1)
        else:
            specie = "Unknown"

        if specie not in out_dict:
            out_dict[specie] = {'ids': set(), 'gc_medio': 0, 'count': 0}

        out_dict[specie]['ids'].add(id)

        comp = seqComposition(seq, "DNA")
        total = sum(comp.values())

        # só conta se for DNA válido (ou seja, só ATCG)
        if total == len(seq):
            gc = (comp['G'] + comp['C']) / total
            out_dict[specie]['gc_medio'] += gc
            out_dict[specie]['count'] += 1

    for specie in out_dict:
        if out_dict[specie]['count'] > 0:
            out_dict[specie]['gc_medio'] = round(
                out_dict[specie]['gc_medio'] / out_dict[specie]['count'], 4
            )
        else:
            out_dict[specie]['gc_medio'] = None

        del out_dict[specie]['count']

    return out_dict


# Tests:

#T1
# print(f"CY043485.1.fasta\n{readFASTA('data_sequences/CY043485.1.fasta')}")
# d =readFASTA("data_sequences/U49845.1.fasta")
# print(d)
# for k, v in d.items(): # v- tuple # k-str
#     print(type(k), type(v), len(v))

# from Bio import SeqIO # problema: ModuleNotFoundError: No module named 'Bio'    solucao:  python -m pip install biopython   OU  pip install biopython 

# bio_dict = {}
# for record in SeqIO.parse("data_sequences/U49845.1.fasta", "fasta"):
#     bio_dict[record.id] = (record.description[len(record.id):].strip(), str(record.seq))

# my_dict = readFASTA("data_sequences/U49845.1.fasta")

# print(my_dict == bio_dict)

#T2
# print(seqComposition('AGCCGTACGG', 'DNA'))
# print(Counter('AGCCGTACGG')) # comparação de contagem
# print(seqComposition('AGCCGTACGG', 'RNA'))
# print(seqComposition('AGCXCGTA_CGGN', 'DNA')) # tem de dar igual ao primeiro teste
# print(seqComposition('AGCXCG_NTACGG', 'NO')) # wriong type
# print(seqComposition('ARNDCEQHILKMFPSTWYV', 'PROTEIN'))

#T2
# print(dnaGCcontent("AGCCGTACGG"))
# print(dnaGCcontent("A"))
# print(dnaGCcontent("G"))
# print(dnaGCcontent("AGCCGXXTACGG")) # uso de SeqComposition caracteres maus, nao afetam GC content
# from Bio.SeqUtils import gc_fraction
# tool_gc = gc_fraction("AGCCGTACGG")
# print(round(tool_gc * 100,2))

#T3
# dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
# rna_seq = transcribeDNA2RNA(dna_seq)
# print("DNA:", dna_seq)
# print("RNA:", rna_seq)
# dna_seq = "AAGGCC"
# rna_seq = transcribeDNA2RNA(dna_seq)
# print("DNA2:", dna_seq)
# print("RNA2:", rna_seq)
# dna_seq = "TTTTT"
# rna_seq = transcribeDNA2RNA(dna_seq)
# print("DNA3:", dna_seq)
# print("RNA3:", rna_seq)
# dna_seq = ""
# rna_seq = transcribeDNA2RNA(dna_seq)
# print("DNA4:", dna_seq)
# print("RNA4:", rna_seq)
# dna_seq = "AAG_TGCTXC" # nada diz sobre isto no enunciado portanto aqui vai so ignorar e retorna-los
# rna_seq = transcribeDNA2RNA(dna_seq)
# print("DNA5:", dna_seq)
# print("RNA5:", rna_seq)

#T4
# nao parar quando encotra codao STOP

# rna_seq = "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG"
# protein_seq = translateRNA2protein(rna_seq)
# print("RNA:", rna_seq)
# print("Protein:", protein_seq)
# rna_seq = "AAGXXXXAUG"
# # nada sobre bases erradas, so ignora mas aqui vejo como ignora tambem nao ter string multipla de tres (sem dar erro): ignora bases a mais no fim
# protein_seq= translateRNA2protein(rna_seq)
# print("RNA2:", rna_seq)
# print("Protein2:", protein_seq)
# rna_seq = ""
# protein_seq = translateRNA2protein(rna_seq)
# print("RNA3:", rna_seq)
# print("Protein3:", protein_seq)

# from Bio.Seq import Seq

# rna_seq = "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG"
# print(Seq(rna_seq).translate())

#T5
# dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCGATAG"
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement:", rev_complement)
# from Bio.Seq import Seq
# print(f"Reverse complement: {str(Seq(dna_seq).reverse_complement())}")
# dna_seq = "ATcGAT"
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement2:", rev_complement)
# dna_seq = "A"
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement3:", rev_complement)
# dna_seq = "CG" # dar CG
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement4:", rev_complement)
# dna_seq = "" # dar CG
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement5:", rev_complement)
# dna_seq = "CXXG" # dar CG; ignora bases más devido ao uso direto do dict
# rev_complement = reverseComplement(dna_seq)
# print("Reverse complement6:", rev_complement)

#T6
# dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
# motif = "GCC"
# positions = findMotif(dna_seq, motif)
# print(positions)
# print(findMotif("GCCATG", "GCC"))
# print(findMotif("ATGGCC", "GCC"))
# print(findMotif("GCAATG", "GCC")) # nao existe
# print(findMotif("AAAAA", "AAA")) # overlap
# print(findMotif("GCCA", "GCCA"))
# print(findMotif("GCC", "GCCC")) #empty list


#T7
# dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
# k = 3
# kmers = mostFrequentKMotifs(dna_seq, k)
# print(kmers)
# dna_seq2 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCGATAG"
# kmers2 = mostFrequentKMotifs(dna_seq2, k)
# print(kmers2)
# print(mostFrequentKMotifs("AAAAAA", 2))
# print(mostFrequentKMotifs("ATATAT", 2))
# print(mostFrequentKMotifs("AAGGTT", 1))
# print(mostFrequentKMotifs("AAGGTT", 1))
# print(mostFrequentKMotifs("ATG", 5))
# print(mostFrequentKMotifs("", 3))
# print(mostFrequentKMotifs("ATG", 3))

#T8
#print(highestGC("data_sequences/J02459.1.fasta"))
#print(highestGC("data_sequences/P68871.2.fasta"))

#T9
# print(compositionMatrix("data_sequences/U49845.1.fasta"))

# result = compositionMatrix("data_sequences/U49845.1.fasta")
# print(compositionMatrix("test.fasta"))

#T10
# dict1 = indexSpeciesGC(readFASTA("data_sequences/NP_001138820.1.fasta"))
# print(dict1)
# dict2 = indexSpeciesGC(readFASTA("data_sequences/NP_001362750.1.fasta"))
# print(dict2)

print("\nT10 structure test:")

d = indexSpeciesGC(readFASTA("data_sequences/NP_001138820.1.fasta"))
print(d)
for specie, info in d.items():
    print(type(specie), type(info), info.keys())

print("\nT10 protein test:")

d = indexSpeciesGC(readFASTA("data_sequences/P68871.2.fasta"))
print(d)
