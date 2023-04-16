#%% 1.2.6
def Text(t, i, k):
    return(t[i:i+k])
def PatternCount(t, Pattern):
    count = 0
    for i in range(len(t) - len(Pattern) + 1):
        Test = Text(t, i, len(Pattern))
        if Test == Pattern:
            count += 1
    return(count)

f = open('dataset_925783_6 (5).txt')
txt = f.readline().strip()
ptn = f.readline().strip()

PatternCount(txt, ptn)

#%% 1.3.2
def ReverseComplement(s):
    bases = {
        "A" : "T",
        "T": "A",
        "C" : "G",
        "G" : "C"        
        }
    rev = ""
    for i in range(len(s)):
        c = s[i]
        rev = bases[c] + rev
    return(rev)

f = open('dataset_925784_2.txt')
txt = f.readline().strip()
ReverseComplement(txt)

def PatternMatch(t, Pattern):
    indices = []
    for i in range(len(t) - len(Pattern) + 1):
        Test = Text(t, i, len(Pattern))
        if Test == Pattern:
            indices.append(str(i))
    print(len(indices), "matches")
    result = ", ".join(indices)
    return(result)

f = open('Vibrio_cholerae.txt')
txt = f.readline().strip()
PatternMatch(txt, "ATGATCAAG")

#%% 1.4.5
def FrequencyTable(t, k):
    freqMap = {}
    for i in range(len(t) - k + 1):
        Test = Text(t, i, k)
        if Test in freqMap:
            freqMap[Test] += 1
        else:
            freqMap[Test] = 1
    return(freqMap)
def FindClumps(t, k, L, T):
    Patterns = []
    n = len(t)
    for i in range(n - L + 1):
        window = Text(t, i, L)
        freqMap = FrequencyTable(window, k)
        for key in freqMap.keys():
            if freqMap[key] >= T:
                Patterns.append(key)
    result = list(set(Patterns)) # remove duplicates
    return(" ".join(result))


f = open('dataset_925785_5 (1).txt')
txt = f.readline().strip()
inputs = f.readline().strip()
FindClumps(txt, 10, 103, 4)

#%% 1.7.10
def Skew(t):
    totalC = 0
    totalG = 0
    diffs = []
    for i in range(len(t)):
        if t[i] == "C":
            totalC += 1
        elif t[i] == "G":
            totalG += 1
        diff = totalG - totalC
        diffs.append(diff)
    minDiff = min(diffs)
    indices = []
    for i in range(len(diffs)):
        if diffs[i] == minDiff:
            indices.append(str(i + 1))
    return(" ".join(indices))
            
f = open('dataset_925788_10(1).txt')
txt = f.readline().strip()
Skew(txt)

#%% 1.8.3
def Hamming(t, u):
    h = 0
    for i in range(len(t)):
        if t[i] != u[i]:
            h += 1
    return(h)

f = open('dataset_925789_3.txt')
txt1 = f.readline().strip()
txt2 = f.readline().strip()
Hamming(txt1, txt2)

#%% 1.8.6
def ApproxMatch(t, Pattern, d):
    indices = []
    for i in range(len(t) - len(Pattern) + 1):
        Test = Text(t, i, len(Pattern))
        if Hamming(Test, Pattern) <= d:
            indices.append(str(i))
    print(len(indices), "matches")
    result = ", ".join(indices)
    return(result)
def ApproxCount(t, Pattern, d):
    count = 0
    for i in range(len(t) - len(Pattern) + 1):
        Test = Text(t, i, len(Pattern))
        if Hamming(Test, Pattern) <= d:
            count += 1
    return(count)


f = open('dataset_925789_6.txt')
pattern = f.readline().strip()
txt = f.readline().strip()
d = int(f.readline().strip())
ApproxCount(txt, pattern, d)



#%% 1.8.9
def FrequentWordsWithMismatches(t, k, d):
    freqMap = {}
    patterns = []
    for i in range(len(t) - k + 1):
        Test = Text(t, i, k)
        Neighborhood = Neighbors(Test, d)
        for n in Neighborhood:
            if n in freqMap:
                freqMap[n] += 1
            else:
                freqMap[n] = 1
    maxMatch = max(freqMap.values())
    for key in freqMap.keys():
        if freqMap[key] == maxMatch:
            patterns.append(key)
    return(patterns)

f = open('dataset_925789_9.txt')
txt = f.readline().strip()
inputs = f.readline().strip()
k = int(inputs[0])
d = int(inputs[2])
FrequentWordsWithMismatches(txt, k, d)
    

#%% 1.11.4
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

f = open('dataset_925792_4.txt')
txt = f.readline().strip()
d = int(f.readline().strip())
" ".join(Neighbors(txt, d))
