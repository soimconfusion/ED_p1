## TODO: Include here the Python standard library modules you need.
#import Biopython #just to use for checks

## Part I

### T1

def readFASTA(file: str):
    with open(file, 'r',encoding='utf-8-sig') as f:
        d = {} # da para pre definir como vai ser? nao, fazer so depois
        ID = None
        seq = ""
        desc= ""
        for line in f:
            line = line.rstrip("\n")
            if line.startswith('>'): # header
                if ID is not None: # ja lemos sequencia daqui
                    d[ID] = (desc, seq)
                #else: 1 header ou qualquer novo
                #ate primeiro espaço
                #index = line.find(" ") nao usar porque caso tenha so a parte do idseq sem descriçao, nao ha espaço da -1
                header = line[1:].strip()
                spl= header.split(maxsplit=1) #primeiro espaço
                ID = spl[0]
                if len(spl) > 1: desc = spl[1]
                else: desc=""
                seq = ""
            else: # sequence
                seq += line # ja tirou-se quebra de linha
                # fasta files nao tem espaços e outras coisas
            #parse last seq
            if ID is not None: # ja lemos sequencia daqui
                    d[ID] = (desc, seq)
    return d

#print(readFASTA('data_sequences/CY043485.1.fasta'))
#print(readFASTA('test.fasta'))

### T2

def seqComposition(seq, type):
    if type == "DNA":
        dict = {'A':0, 'T':0, 'C':0, 'G':0}
    
    elif type == "RNA":
        dict = {'A':0, 'U':0, 'C':0, 'G':0}
    
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
#print(seqComposition('AGCCGTACGG', 'DNA'))
#print(seqComposition('AGCCGTACGG', 'RNA'))

def dnaGCcontent(seq):
   dict = seqComposition(seq, "DNA")
   cg = dict['C'] + dict['G']
   total = len(seq)
   return str(round(cg/total * 100)) + ('%'+ " GC composition")

#print(dnaGCcontent("AGCCGTACGG"))
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

dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
rna_seq = transcribeDNA2RNA(dna_seq)
print("DNA:", dna_seq)
print("RNA:", rna_seq)
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
            genetic_code.append(temp2) # AAs; start; base_1; base_2; base3;
        AAs = genetic_code[0]
        base_1 = genetic_code[2]
        base_2 = genetic_code[3]
        base_3 = genetic_code[4]
        for i in range(0, len(base_1)-1):
            dict[base_1[i]+base_2[i]+base_3[i]] = AAs[i]
        for j in range(0, len(seq), 3):
            codon = str(seq[j] + seq[j+1] + seq[j+2])
            if codon in dict.keys():
                protein += dict[codon]
        return protein
    
rna_seq = "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG"
protein_seq = translateRNA2protein(rna_seq)
print("RNA:", rna_seq)
print("Protein:", protein_seq)


### T5

def reverseComplement(sequence : str):
    #complement:
    dna_complement = {'A': 'T', 'C':'G', 'T': 'A', 'G':'C'}
    dna_com = []
    for base in sequence: #0(n)
        if base in dna_complement.keys():
            dna_com.append(dna_complement[base])
    return "".join(dna_com[::-1]) # or use list.reverse(); reversed()<-iterator or use a loop
                #join -> O(n)
dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCGATAG"
rev_complement = reverseComplement(dna_seq)
print("Reverse complement:", rev_complement)

## Parte III

### T6

def findMotif(dna_sequence : str, motif :str):
    step = len(motif)
    positions = []
    find = []
    # most be a better way to do this...
    for i in range(0, len(dna_sequence) - step + 1):
        find.append(dna_sequence[i])
        if len(find) == 3:
            if "".join(find) == motif:
                positions.append(i - 2)
            find = find[1:]
    return positions

dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
motif = "GCC"
positions = findMotif(dna_seq, motif)
#print(positions)

### T7
from collections import Counter
def mostFrequentKMotifs(dna_sequence: str, k : int):
    motifs = []
    result = []
    for b in range(0, len(dna_sequence) -k +1): #o(n)
            motifs.append(dna_sequence[b: b + k]) #O(k) ?
    most_motifs = Counter(motifs)
    max_freq = most_motifs.most_common(1)[0][1]
    for i,y in most_motifs.items():
        if y == max_freq:
            result.append((i,y))
    return result

dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
dna_seq2 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCGATAG"
k = 3
kmers2 = mostFrequentKMotifs(dna_seq2, k)
kmers = mostFrequentKMotifs(dna_seq, k)
#print(kmers)
#print(kmers2)

## Parte IV

### T8

def highestGC(file):
    dict = readFASTA(file)
    ID = list(dict.keys())
    CG_list = []
    seqs = []
    for _, y in dict.values():
        seqs.append(y)
    for seq in seqs:
        type_count = Counter(seq)
        cg = type_count['C'] + type_count['G']
        total = len(seq)
        CG_list.append(round(cg/total * 100))
    i = CG_list.index(max(CG_list))
    return (ID[i], CG_list[i])
#highestGC("data_sequences/J02459.1.fasta")
#print(highestGC("data_sequences/J02459.1.fasta"))
### T9

#trocar para DNA
def compositionMatrix(file):
    dict = readFASTA(file)
    ID = list(dict.keys())
    dict_interno = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    lis =[]
    result={}
    for _, y in dict.values():
        for i in range(len(y)):
            if y[i] in dict_interno.keys():
                dict_interno[y[i]] += 1
        lis.append(dict_interno)
    
    for l in range(0, len(ID)):
        result[ID[l]] = lis[l]
    return result
#print(compositionMatrix("test.fasta"))

### T10
# considerar "proteins" na verdade se for RNA tambem vai tratar so considerar ATCGN!
import re
def indexSpeciesGC (fasta: dict):
    out_dict = {} #tipo dict[][] [species][(ids: set(identifiers), 'gc_medio)]
    # nao é preciso uma lista para cada passa/info
    #recolher toda a info necessaria de fasta de uma so vez:
    
    # falta por o id certo é preciso que seja um set antes de se começara indexar...

    for id, (description, seq) in fasta.items():
        found = re.search(r"\[(.*?)\]", description)
        if found:
            specie = found.group(1)
        else:
            specie = "Unknown"
        if specie not in out_dict: #new key
            out_dict[specie] = {'ids': set(), 'gc_medio': 0, 'count': 0}
        out_dict[specie]['ids'].add(id) # ja sabe que é set()
# gc_count é None caso tenha algum caracter alem de `A`, `T`, `C`, `G` or `N`
        valid_seq = True
        #count = 0
        for c in seq:
            if c not in "ATCGN":
                valid_seq = False
        if valid_seq:
            #na verdade toal aqui!
            out_dict[specie]['gc_medio'] += (seq.count('G') + seq.count('C')) / len(seq)
            out_dict[specie]['count'] += 1
            #count += 1 <- store in dict then delete value!
    # ja vimos todas as sequencias!
    for specie in out_dict:
        if out_dict[specie]['count']> 0:
            out_dict[specie]['gc_medio'] = round(out_dict[specie]['gc_medio'] / out_dict[specie]['count'], 4)
        else:
            out_dict[specie]['gc_medio'] = None
        del out_dict[specie]['count']
    
    # DO search if theres a python funct that grabs substrings by chars
        #re.search()<. usa uma sintax tipo sed(?ig)
        # sim usa expressoes regulares (regex)
        # \ para escapar [ e ] como caracteres especiais
    return out_dict
# tem header com o formato bem, mas é proteina
dict1 = indexSpeciesGC(readFASTA("data_sequences/NP_001138820.1.fasta"))
print(dict1)
dict2 = indexSpeciesGC(readFASTA("data_sequences/NP_001362750.1.fasta"))
print(dict2)
#ids_gc_medio  =  { }
# {dict.values[0] : ids_gc_medio }
