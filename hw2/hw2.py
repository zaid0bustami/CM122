#%% 9.3.4
class Node:
    def __init__(self, name = "", edges = []):
        self.name = name
        self.edges = edges
    def setName(self, name):
        self.name = name
    def addEdge(self, edge):
        self.edges.append(edge)
    def getName(self):
        return(self.name)
    def getEdges(self):
        return(self.edges)
    def getEdgesNames(self):
        edgesNames = [i.getName() for i in self.edges]
        return(edgesNames)
class Edge:
    def __init__(self, name = "", start = None, end = None):
        self.name = name
        self.start = start
        self.end = end
        # add edge to the starting node
        # limitation - you can only walk down the edges
        start.addEdge(self)
    def setName(self, name):
        self.name = name
    def getName(self):
        return(self.name)
    def getStart(self):
        return(self.start)
    def getEnd(self):
        return(self.end)
class Trie:
    def __init__(self, Patterns):
        nodeIndex = 0
        self.root = Node(nodeIndex)
        for Pattern in Patterns:
            currentNode = self.root
            for i in range(len(Pattern)):
                currentSymbol = Pattern[i]
                edgesNames = currentNode.getEdgesNames()
                if currentSymbol in edgesNames:
                    e = edgesNames.index(currentSymbol)
                    edge = currentNode.getEdges()[e]
                    currentNode = edge.getEnd()
                else:
                    nodeIndex += 1
                    newNode = Node(nodeIndex, [])
                    newEdge = Edge(currentSymbol, currentNode, newNode)
                    currentNode = newNode
                    print(newEdge.getStart().getName(), newEdge.getEnd().getName(), newEdge.getName())
        self.size = nodeIndex + 1
        print("Created trie of size", self.size)
    def getRoot(self):
        return(self.root)

import sys
sys.stdout = open('hw2/output.txt', 'w')

f = open('hw2/dataset_925947_4.txt')
patterns = f.readline().strip().split(sep = " ")
t = Trie(patterns)

sys.stdout.close()

#%% 9.5.4 EXTRA CREDIT



#%% 9.5.5 EXTRA CREDIT



#%% 9.6.2
class SuffixArray:
    def __init__(self, Text):
        Array = []
        for i in range(len(Text)):
            entry = (i, Text[i:len(Text)])
            Array.append(entry)
        self.Array = sorted(Array, key=lambda x: x[1])
        self.Indices = [i[0] for i in self.Array]
    def getIndices(self):
        return(self.Indices)

f = open('hw2/dataset_925950_2.txt')
txt = f.readline().strip()
test = SuffixArray(txt)
' '.join(map(str, test.getIndices()))

#%% 9.7.5
def BWT(Text):
    
    

#%% 9.9.11

    

#%% 9.10.8

    

#%% 9.11.7
    

