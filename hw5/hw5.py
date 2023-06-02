#%% 3.2.3
def Composition(Text, k):
    result = []
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        result.append(kmer)
    return(result)

with open('hw5/dataset_925820_3.txt') as f:
    k = int(f.readline().strip())
    txt = f.readline().strip()
    with open('hw5/result.txt', 'w') as g:
        print(' '.join(Composition(txt, k)), file = g)

#%% 3.3.3
def PathToGenome(path):
    result = path[0]
    for i in range(1, len(path)):
        kmer = path[i]
        result += kmer[-1]
    return(result)

with open("hw5/dataset_925821_3.txt") as f:
    path = f.readline().strip().split(sep = " ")
    print(PathToGenome(path))

#%% 3.3.10
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

def Overlap(Patterns):
    G = Graph()
    for i in range(len(Patterns)):
        for j in range(len(Patterns)):
            Pattern1 = Patterns[i]
            Pattern2 = Patterns[j]
            if Pattern1[1:] == Pattern2[:-1]:
                G.newEdge(Pattern1, Pattern2)
    return(G)

with open('hw5/dataset_925821_10.txt') as f:
    patterns = f.readline().strip().split(sep = " ")
    with open('hw5/result.txt', 'w') as g:
        print(Overlap(patterns).writeEdges(file=g))

#%% 3.4.6
def DeBruijn(Text, k):
    kmers = Composition(Text, k - 1)
    result = Graph()
    for i in range(len(kmers) - 1):
        result.newEdge(kmers[i], kmers[i+1])
        result.newEdge
    return(result)

with open('hw5/dataset_925822_6.txt') as f:
    k = int(f.readline().strip())
    txt = f.readline().strip()
    with open('hw5/output.txt', 'w') as g:
        DeBruijn(txt, k).writeEdges(file = g)
        test = DeBruijn(txt, k)

#%% 3.5.8
def DeBruijn2(Patterns):
    result = Graph()
    for i in range(len(Patterns)):
        kmer = Patterns[i]
        result.newEdge(kmer[:-1], kmer[1:])
    return(result)

with open('hw5/dataset_925823_8.txt') as f:
    patterns = f.readline().strip().split(' ')
    with open('hw5/output.txt', 'w') as g:
        DeBruijn2(patterns).writeEdges(file = g)

#%% 3.8.2
import random

def TextToAdjList(Text):
    result = dict()
    for i in Text:
        i = i.split(sep=": ")
        j = {int(i[0]) : i[1]}
        j[int(i[0])] = [int(k) for k in j[int(i[0])].split(" ")]
        result.update(j)
    return(result)

def AdjListToGraph(AdjList):
    G = Graph()
    G.setRoot(list(AdjList.keys())[0])
    for key in AdjList:
        for value in AdjList[key]:
            G.newEdge(key, value)
    return(G)

# def EulerianCycle(Graph):
#     EdgesLeft = Graph.Edges.copy()
    
#     # generate all paths
#     paths = []
#     while len(EdgesLeft) > 0:  # change later
#         # the first step is unique bc of the randomly chosen node
#         path = []
#         Nodes = list(EdgesLeft.keys())
#         prevNode = random.choice(Nodes)
#         path.append(prevNode)
#         nextNode = EdgesLeft[prevNode][0]
#         EdgesLeft[prevNode].remove(nextNode)
#         if EdgesLeft[prevNode] == []:
#             del EdgesLeft[prevNode]
#         # the rest of the steps are the same
#         while nextNode in EdgesLeft:
#             prevNode = nextNode
#             path.append(prevNode)
#             nextNode = EdgesLeft[prevNode][0]
#             EdgesLeft[prevNode].remove(nextNode)
#             if EdgesLeft[prevNode] == []:
#                 del EdgesLeft[prevNode]
#         paths.append(path)
        
#     # make the first path a cycle and connect the paths
#     paths[0].append(paths[0][0])
#     print(paths)
#     for p in reversed(range(len(paths) - 1)):
#         firstPath = paths[p]
#         secondPath = paths[p + 1]
#         commonNode = set(firstPath).intersection(set(secondPath)).pop()
        
#         firstIndex = firstPath.index(commonNode) + 1
#         secondIndex = secondPath.index(commonNode) + 1
        
#         if firstIndex > len(firstPath) - 1:
#             firstPath.append(None)
#         if secondIndex > len(secondPath) - 1:
#             secondPath.append(None)
#         firstNext = firstPath[firstIndex]
#         secondNext = secondPath[secondIndex]
        
#         firstPath[firstIndex] = secondNext
#         secondPath[secondIndex] = firstNext
    
#     # create the cycle
#     cycle = []
#     for p in paths:
#         for n in p:
#             if n != None:
#                 cycle.append(n)
#     return(cycle)

def EulerianCycle(Graph):
    cycle = []
    start_node = list(Graph.keys())[0]  # Start from any node in the graph
    print(start_node)
    def explore_path(node):
        while node in Graph and Graph[node]:  # While there are remaining edges from the current node
            next_node = Graph[node].pop(0)  # Remove the first edge from the current node
            explore_path(next_node)  # Recursively explore the next node
        cycle.append(node)  # Add the current node to the cycle

    explore_path(start_node)  # Start exploring from the start node

    cycle.reverse()  # Reverse the cycle to obtain the Eulerian cycle
    print(Graph)
    return cycle


with open('hw5/dataset_925826_2.txt') as f:
    Text = f.readlines()
    AdjList = TextToAdjList(Text)
    G = AdjListToGraph(AdjList)
    cycle = EulerianCycle(G.Edges)
    print(" ".join([str(i) for i in cycle]))

#%% 3.8.6
# did by hand

#%% 3.8.7
import networkx as nx

def form_string_from_kmers(k, kmers):
    G = nx.MultiDiGraph()

    # Add edges to the graph
    for kmer in kmers:
        # Nodes are (k-1)mers, edges are kmers
        G.add_edge(kmer[:k-1], kmer[1:], kmer=kmer)

    # Find the Eulerian path
    path = list(nx.eulerian_path(G))

    # Start with the first node of the first edge
    result = path[0][0]

    # Then add the end node of each edge
    for _, end_node in path:
        result += end_node[-1]

    return result

with open('hw5/dataset_925826_7.txt') as f:
    k = int(f.readline().strip())
    stringList = f.readline().strip().split(" ")
    print(form_string_from_kmers(k, stringList))

