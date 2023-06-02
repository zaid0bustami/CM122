#%% Libraries
import os
import statistics as st

#%% Functions
def readsFASTA(path):
    with open(path) as f:
        reads = []
        names = []
        currentRead = []
        for line in f.readlines():
            if line[0] == ">":
                if len(currentRead) != 0:
                    reads.append(''.join(currentRead))
                    currentRead = []
                names.append(line.strip())
            else:
                currentRead.append(line.strip())
        result = []
        for i in range(len(reads)):
            result.append((names[i], reads[i]))
        return(result)

def genomeFASTA(path):
    with open(path) as f:
        meta = f.readline().strip().split(" ")
        name = meta[0]
        length = meta[2]
        genome = ''.join([i.strip() for i in f.readlines()])
        result = (name, length, genome)
        return(result)
    
def HashTable(genome, k):
    result = dict()
    for i in range(len(genome) - k + 1):
        kmer = genome[i : i + k]
        if kmer in result:
            result[kmer].append(i)
        else:
            result[kmer] = [i]
    return(result)

def MatchRead(index, read, k):
    result = dict()
    result[read] = []
    if len(read) < k:
        result[read].append([-1])
    else:
        for i in range(len(read) - k + 1):
            kmer = read[i : i + k]
            if kmer in index:
                pos = index[kmer]
            else:
                pos = [-1]
            result[read].append(pos)
    return(result)

def MapReads(reads, genomeHashTable, k):
    readMap = dict()
    for i in range(len(reads)):
        read = reads[i][1]
        match = MatchRead(genomeHashTable, read, k)
        allMatches = [j - i for i in range(len(match[read])) for j in match[read][i] if j != -1]
        if len(allMatches) == 0:
            continue
        result = {read: st.mode(allMatches)}
        readMap.update(result)    
    return(readMap)

def MapGenomes(genomes, reads, k):
    print("Mapping to Genomes")
    result = []
    for i in range(len(genomes)):
        if i % 10 == 0:
            print(".", end="")
        name = genomes[i][0]
        genome = genomes[i][2]
        genomeHashTable = HashTable(genome, k)
        readMap = MapReads(reads, genomeHashTable, k)
        add = (name, len(readMap))
        result.append(add)
    print(".")
    return(result)

#%% 4a
#%% Inputs
path = "data\\project4a-data\\"
reads = readsFASTA(path + "project4a_10000_reads.fasta")
genomeFiles = [path + f for f in os.listdir(path) if "genome" in f]
genomes = [genomeFASTA(f) for f in genomeFiles]

#%% Read-Mapping for Each Genome
k = 31
genomeMap = sorted(MapGenomes(genomes, reads, k), key = lambda x : x[1], reverse=True)
matchedGenomes = [i[0] for i in genomeMap if i[1] > 500]

#%% Output
with open("predictions.csv", 'w') as f:
    readMaps = []
    gs = [g[0] for g in genomes]
    for g in matchedGenomes:
        print(".", end = "")
        index = gs.index(g)
        genome = genomes[index][2]
        genomeHashTable = HashTable(genome, k)
        readMap = MapReads(reads, genomeHashTable, k)
        readMaps.append((g, readMap))
    for read in reads:
        matched = False
        for readMap in readMaps:
            if read[1] in readMap[1]:
                matched = True
                print(read[0] + "\t" + readMap[0][1:], file = f)
        if matched == False:
            print(read[0] + "\t" + matchedGenomes[0][1:], file = f)
    print('>read_19999' + "\t" + matchedGenomes[0][1:], file = f)

                
            
        