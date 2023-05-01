#%% 2.2.8
def Text(t, i, k):
    return(t[i:i+k])

def Hamming(t, u):
    h = 0
    for i in range(len(t)):
        if t[i] != u[i]:
            h += 1
    return(h)

def Neighbors(t, d):
    bases = {"A", "C", "G", "T"}
    if d == 0:
        return({t})
    if len(t) == 1:
        return(bases)
    neighborhood = set()
    suffixT = t[1:len(t)]
    suffixNeighbors = Neighbors(suffixT, d)
    for s in suffixNeighbors:
        if Hamming(suffixT, s) < d:
            for b in bases:
                neighborhood.add(b + s)
        else:
            neighborhood.add(t[0] + s)
    return(neighborhood)

def ApproxCount(t, Pattern, d):
    count = 0
    for i in range(len(t) - len(Pattern) + 1):
        Test = Text(t, i, len(Pattern))
        if Hamming(Test, Pattern) <= d:
            count += 1
    return(count)

def MotifEnumeration(Dna, k, d):
    Patterns = set()
    for Sequence in Dna:
        print(Sequence)
        for i in range(len(Sequence) - k + 1):
            Pattern = Sequence[i:i+k]
            for Pattern2 in Neighbors(Pattern, d):
                test = [ApproxCount(j, Pattern2, d) > 0 for j in Dna] # test if there is an approx match for Pattern2 for every sequence in Dna
                if sum(test) / len(test) == 1: # if every value in test is true
                    Patterns.add(Pattern2)
    return(Patterns)

with open('hw4/dataset_925801_8.txt') as f:
    inputs = [int(i) for i in f.readline().strip().split(' ')]
    txt = f.readline().strip().split(' ')
' '.join(list(MotifEnumeration(txt, inputs[0], inputs[1])))

#%% 2.5.2
def ProfileMostProb(Text, k, Profile):
    maxProb = (0, Text[0 : k])
    index = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    for i in range(len(Text) - k + 1):
        kmer = Text[i : i + k]
        prob = 1
        for j in range(k):
            c = kmer[j]
            row = index[c]
            prob *= Profile[row][j]
        if prob > maxProb[0]:
            maxProb = (prob, kmer)
    return(maxProb[1])
            
with open('hw4/dataset_925804_3.txt') as f:
    txt = f.readline().strip()
    k = int(f.readline().strip())
    mat = []
    for row in f.readlines():
        entry = [float(i) for i in row.strip().split(sep=' ')]
        mat.append(entry)

ProfileMostProb(txt, k, mat)
                
#%% 2.5.3
0

#%% 2.5.5
import numpy as np

def Profile(Motifs, k):
    index = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    counts = np.zeros((4,k))
    for Motif in Motifs:
        for i in range(k):
            c = Motif[i]
            row = index[c]
            counts[row, i] += 1
    n = len(Motifs)
    counts[counts != 0] /= n
    return(counts)

def Score(Motifs, Profile, k):
    score = 0
    optimal = ""
    index = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
    for i in range(k):
        probs = [row[i] for row in Profile]
        maxProb = max(probs)
        optimalBase = index[probs.index(maxProb)]
        optimal += optimalBase
        actualBases = [row[i] for row in Motifs]
        for actualBase in actualBases:
            if actualBase != optimalBase:
                score += 1                
    return(score)
    
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = [Seq[0 : k] for Seq in Dna]
    for i in range(len(Dna[0]) - k + 1):
        Motif = Dna[0][i : i + k]
        Motifs = [Motif]
        for j in range(1, t):
            profile = Profile(Motifs[0 : j], k)
            Motifs.append(ProfileMostProb(Dna[j], k, profile))
        if Score(Motifs, Profile(Motifs, k), k) < Score(BestMotifs, Profile(BestMotifs, k), k):
            BestMotifs = Motifs
    return(BestMotifs)

with open('hw4/dataset_925804_5.txt') as f:
    inputs = [int(i) for i in f.readline().strip().split(' ')]
    dna = f.readline().strip().split(sep=' ')

" ".join(GreedyMotifSearch(dna, inputs[0], inputs[1]))


#%% 2.7.5 Extra Credit


#%% 2.9.4 Extra Credit