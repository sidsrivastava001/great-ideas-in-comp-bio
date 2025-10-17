import random

"""
Motif finding using Gibbs sampling and random k-mer search
"""

"""
Create a k-mer frequency table using length k windows
Takes in DNA string, k-mer length k
Returns dict with k-mers as keys, frequency as values
"""
def freqTable(DNA, k):
    d = dict()
    for i in range(len(DNA)-k+1):
        if (DNA[i:i+k] not in d):
            d[DNA[i:i+k]] = 1
        else:
            d[DNA[i:i+k]] += 1
    return d


"""
Look for most probable k-mer in DNA string given nucleotide frequency table
Takes in DNA string, k-mer length k
Returns most probable k-mer
"""
def mostProbableKmer(DNA, k, matrix):
    prob = 0
    most = ""
    for i in range(len(DNA)-k):
        p = 1
        for j in range(i, i+k):
            p*=float(matrix[DNA[j]][j-i])
        if(p>=prob):
            most = DNA[i:i+k]
            prob = p
    return most

"""
Gets score of list of motifs
Takes in list of Motifs
Returns score value
"""
def score(Motifs):
    score = 0
    max = 0
    maxchar = ""
    d = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for i in range(len(Motifs[0])):
        for j in range(len(Motifs)):
            d[Motifs[j][i]]+=1
            if(d[Motifs[j][i]]>max):
                max = d[Motifs[j][i]]
                maxchar = Motifs[j][i]
        d = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for j in range(len(Motifs)):
            if(Motifs[j][i] != maxchar):
                score+=1
        max = 0
        maxchar = ""
    return score

"""
Create a profile matrix for a set of motifs based on nucleotide frequency
Takes in list of Motifs
Returns Profile matrix with size 4 x motif size
"""
def Profile(Motifs):
    d = {'A': [], 'T': [], 'C': [], 'G': []}
    for i in range(len(Motifs[0])):
        d['A'].append(1)
        d['T'].append(1)
        d['C'].append(1)
        d['G'].append(1)
    for i in range(len(Motifs[0])):
        for j in range(len(Motifs)):
            d[Motifs[j][i]][i]+=1
    for i in range(len(Motifs[0])):
        total = d['A'][i] + d['T'][i] + d['C'][i] + d['G'][i]
        d['A'][i] = d['A'][i]/total
        d['T'][i] = d['T'][i]/total
        d['C'][i] = d['C'][i]/total
        d['G'][i] = d['G'][i]/total
    return d

"""
Look for most probably k-mers in DNA string with randomized motif search
Takes in DNA string, k-mer size
Return list of best motifs
"""
def RandomizedMotifSearch(DNA, k):
    BestMotifs = []
    for i in range(len(DNA)):
        x = random.randint(0, len(DNA[i])-k)
        BestMotifs.append(DNA[i][x:x+k])
    Motifs = BestMotifs
    while(1):
        profile = Profile(Motifs)
        Motifs = []
        for i in DNA:
            a = mostProbableKmer(i, k, profile)
            Motifs.append(a)
        if(score(Motifs) < score(BestMotifs)):
            BestMotifs = Motifs
        else:
            return BestMotifs

"""
Generate random k-mer weighted by profile matrix of nucleotides
Takes in DNA string, k-mer length, and profile matrix
Returns random k-mer
"""
def RandomKmer(DNA, k, matrix):
    l = []
    probs = []
    for i in range(len(DNA)-k):
        p = 1
        for j in range(i, i+k):
            p*=float(matrix[DNA[j]][j-i])
        l.append(DNA[i:i+k])
        probs.append(p)
    return random.choices(l, weights=probs)

"""
Gibbs Sampler to find most common motif in DNA string
Takes in DNA string, k-mer length, and n iterations of sampler to run
Returns list of BestMotifs
"""
def GibbsSampler(DNA, k, n):
    BestMotifs = []
    for i in range(len(DNA)):
        x = random.randint(0, len(DNA[i])-k)
        BestMotifs.append(DNA[i][x:x+k])
    Motifs = BestMotifs
    for j in range(n):
        rand = random.randint(0, len(DNA)-1)
        Motifs.pop(rand)
        profile = Profile(Motifs)
        a = RandomKmer(DNA[rand], k, profile)[0]
        Motifs.insert(rand, a)
        if(score(Motifs) < score(BestMotifs)):
            BestMotifs = Motifs
    return BestMotifs