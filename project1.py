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

def reverseComplement():
    return None


## Parte III

### T6

def findMotif():
    return None

### T7

def mostFrequentKMotifs():
    return None


## Parte IV

### T8

def highestGC():
    return None

### T9

def compositionMatrix():
    return None

### T10

def indexSpeciesGC ():
    return None