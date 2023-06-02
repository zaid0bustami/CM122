#%% Libraries
import matplotlib.pyplot as plt
import random
import statistics as st

#%% Classes and Functions  
# Assembly
class Graph:
    def __init__(self, root = None):
        self.Nodes = dict()
        self.Edges = dict()
        self.Root = root
    def setRoot(self, root):
        self.Root = root
    def getRoot(self):
        return(self.Root)
    def newNode(self, one, data = None):
        if one not in self.Nodes:
            self.Nodes[one] = data
    def newEdge(self, one, two, oneData = None, twoData = None):
        self.newNode(one, oneData)
        self.newNode(two, twoData)
        if one in self.Edges:
            self.Edges[one].append(two)
        else:
            self.Edges[one] = [two]
    def displayNodes(self):
        print("Nodes")
        print("ROOT* = " + str(self.Root))
        for i in self.Nodes:
            print(str(i) + " = " + str(self.Nodes[i]))
    def writeEdges(self, file):
        for i in self.Edges:
            result = [str(j) for j in self.Edges[i]]
            print(str(i) + ": " + " ".join(result), file = file)
    def displayEdges(self):
        print("Edges")
        for i in self.Edges:
            result = [str(j) for j in self.Edges[i]]
            print(str(i) + ": " + " ".join(result))
    def display(self):
        print("--------------------")
        self.displayNodes()
        print("--------------------")
        self.displayEdges()

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

def Composition(Text, k):
    result = []
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        result.append(kmer)
    return(result)

def ReadHashTable(reads, k):
    result = dict()
    for i in range(len(reads)):
        read = reads[i][1]
        if i % 10000 == 0:
            print(".", end = "")
        composition = Composition(read, k)
        for kmer in composition:
            if kmer in result:
                result[kmer] += 1
            else:
                result[kmer] = 1
    print(".")
    return(result)

def Spectrum(counts, c1, c2):
    result = []
    for kmer in counts:
        count = counts[kmer]
        if count <= c1:
            repeat = 0
        elif count <= c2:
            repeat = 1
        else:
            repeat = 2
        for j in range(repeat):
            result.append(kmer)
    return(result)

def SpectrumDeBruijn(spectrum):
    result = Graph()
    for kmer in spectrum:
        result.newEdge(kmer[:-1], kmer[1:])
    return(result)

# def InDegree(node, dB):
#     inDegree = 0
#     for i in dB.Edges:
#         if node in dB.Edges[i]:
#             inDegree += 1
#     return(inDegree)

# def EulerDeBruijn(dB): # eulerize the graph by removing odd-degreed nodes
#     result = Graph()
#     result.Edges = dB.Edges.copy()
#     result.Nodes = dB.Nodes.copy()
#     count = 0
#     for i in dB.Edges:
#         outDegree = len(dB.Edges[i])
#         inDegree = InDegree(i, dB)
#         degree = outDegree + inDegree
#         if degree % 2 == 1:
#             if count < 2: # ensure 2 odd-degreed edges
#                 count += 1
#                 continue
#             else:
#                 result.Edges.pop(i)
#                 result.Nodes.pop(i)
#     return(result)

def Paths(dB):
    EdgesLeft = dB.Edges.copy()
    paths = []
    while len(EdgesLeft) > 0:  # change later
        # the first step is unique bc of the randomly chosen node
        path = []
        Nodes = list(EdgesLeft.keys())
        prevNode = random.choice(Nodes)
        path.append(prevNode)
        nextNode = EdgesLeft[prevNode][0]
        EdgesLeft[prevNode].remove(nextNode)
        if EdgesLeft[prevNode] == []:
            del EdgesLeft[prevNode]
        # the rest of the steps are the same
        while nextNode in EdgesLeft:
            prevNode = nextNode
            path.append(prevNode)
            nextNode = EdgesLeft[prevNode][0]
            EdgesLeft[prevNode].remove(nextNode)
            if EdgesLeft[prevNode] == []:
                del EdgesLeft[prevNode]
        paths.append(path)
    return(paths)
    
def Genome(paths):
    result = ""
    for p in paths:
        concatPath = p[0]
        for i in range(len(p) - 1):
            concatPath += p[i + 1][-1]
        result += concatPath
    return(result)

# Read-Mapping
def HashTable(genome, k):
    result = dict()
    for i in range(len(genome) - k + 1):
        if i % 10000 == 0:
            print(".", end = "")
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

#%% 3b
# Analyzing Read Counts
print('Assembly')
reads = readsFASTA("Data\\project3b_20000_reads_without_positions.fasta")
counts = ReadHashTable(reads, k = 20)
c1 = 5
c2 = 39
data = list(counts.values())

fig, axs = plt.subplots(3, 1)
axs[0].set_title("Histogram of Kmer Counts in Reads")
axs[1].set_ylabel("Frequency")
axs[2].set_xlabel("Count")
axs[0].hist([i for i in data if i <= c1],               bins = c1)
axs[1].hist([i for i in data if i > c1 and i <= c2],    bins = c2-c1)
axs[2].hist([i for i in data if i > c2],                bins = max(data)-c2)
plt.show()

# Building Spectrum and Genome
spectrum = Spectrum(counts, c1, c2)
dB = SpectrumDeBruijn(spectrum)
paths = Paths(dB)
genome = Genome(paths)

# Genome Index
print('Storing HashTable')
k = 8
genomeHashTable = HashTable(genome, k)

# Read Mapping
print(".")
print('Mapping Reads')
readMap = []
for i in range(len(reads)):
    label = reads[i][0]
    read = reads[i][1]
    match = MatchRead(genomeHashTable, read, k)
    allMatches = [j - i for i in range(len(match[read])) for j in match[read][i] if j != -1]
    if len(allMatches) == 0:
        continue
    result = (label, read, st.mode(allMatches))
    readMap.append(result)    
    if i % 1000 == 0:
        print(".", end = "")
        
result = [(i[0], i[1]) for i in sorted(readMap, key = lambda x: x[2])]
with open("predictions.csv", 'w') as f:
    for i in result:
        print(i[0], file = f)