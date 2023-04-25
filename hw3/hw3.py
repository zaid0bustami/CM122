#%% Libraries
import math
import numpy as np
from numpy import random as rd

#%% 8.6.2
#calculates the furthest point in DataPoints from DataPoint
def furthestPoint(Centers, DataPoints):
    maxDistance = -1
    furthest = None
    for pt in DataPoints:
        minDistance = float('inf')
        for c in Centers:
            distance = math.dist(c, pt)
            if distance < minDistance:
                minDistance = distance
        if minDistance > maxDistance:
            furthest = pt
            maxDistance = minDistance
    return(furthest)
def FarthestFirstTraversal(k, m, Data):
    i = rd.randint(len(Data))
    Centers = {Data[0]} #replace with i for random initialization
    while len(Centers) < k:
        newCenter = furthestPoint(Centers, Data)
        Centers.add(newCenter)
    return(Centers)

f = open('hw3/dataset_925929_2.txt')
inputs = [int(i) for i in f.readline().strip().split(' ')]
data = []
for line in f.readlines():
    add = [float(i) for i in line.strip().split(' ')]
    data.append(tuple(add))
result = FarthestFirstTraversal(inputs[0], inputs[1], data)
{' '.join(map(str, tup)) for tup in result}

#%% 8.7.3
def Distortion(Data, Centers):
    totalDistance = 0
    for pt in Data:
        closestDistance = float('inf')
        for c in Centers:  
            distance = math.dist(pt, c)
            if distance < closestDistance:
                closestDistance = distance
        totalDistance += closestDistance ** 2
    result = totalDistance / len(Data)
    return(result)

with open('hw3/dataset_925930_3.txt') as f:
    inputs = [int(i) for i in f.readline().strip().split(' ')]
    centers = []
    line = f.readline()
    while line.strip() != '--------':
        add = [float(i) for i in line.strip().split(' ')]
        centers.append(tuple(add))
        line = f.readline()
    data = []
    for line in f.readlines():
        add = [float(i) for i in line.strip().split(' ')]
        data.append(tuple(add))
Distortion(data, centers)

#%% 8.8.3
def closestCenter(Centers, DataPoint):
    closestDistance = float('inf')
    for c in Centers:  
        distance = math.dist(DataPoint, c)
        if distance < closestDistance:
            closestDistance = distance
            closestPoint = c
    return(closestPoint)
def centerOfGravity(m, DataPoints):
    center = [None] * m
    n = len(DataPoints)
    for dimension in range(m):
        total = sum([i[dimension] for i in DataPoints])
        average = total / n
        center[dimension] = average
    return(tuple(center))
def Lloyd(k, m, Data):
    #random initialization
    # Centers = []
    # randints = rd.randint(len(Data), size = k)
    # for i in randints:
    #     Centers.append(Data[i])
    newCenters = Data[0 : k] # remove and uncomment above for random initialization
    Centers = []
    while newCenters != Centers:
        Centers = [i for i in newCenters]
        #assigning datapoints to clusters
        assignments = {}
        for pt in Data:
            closest = closestCenter(Centers, pt)
            assignments[pt] = Centers.index(closest)
        # recalculating centers
        for cluster in range(k):
            clusterPoints = []
            for pt in assignments.keys():
                if assignments[pt] == cluster:
                    clusterPoints.append(pt)
            newCenters[cluster] = centerOfGravity(m, clusterPoints)    
    return(Centers)

with open('hw3/dataset_925931_3.txt') as f:
    inputs = [int(i) for i in f.readline().strip().split(' ')]
    data = []
    for line in f.readlines():
        add = [float(i) for i in line.strip().split(' ')]
        data.append(tuple(add))
Lloyd(inputs[0], inputs[1], data)

#%% 8.14.7
class Graph:
    def __init__(self, root = None):
        self.Nodes = dict()
        self.Edges = dict()
        self.Root = root
    def setRoot(self, root):
        self.Root = root
    def Root(self):
        return(self.Root)
    def newNode(self, one, data = None):
        if one not in self.Nodes:
            self.Nodes[one] = data
            self.Edges[one] = []
    def newEdge(self, one, two):
        self.newNode(one)
        self.newNode(two)
        self.Edges[one].append(two)
    def displayNodes(self):
        print("Nodes")
        for i in self.Nodes:
            print(str(i) + " = " + str(self.Nodes[i]))
    def displayEdges(self):
        print("Edges")
        for i in self.Edges:
            result = [str(j) for j in self.Edges[i]]
            print(str(i) + " --> " + " ".join(result))
def avgDist(D, ones, twos):
    total = 0
    for i in ones:
        for j in twos:
           total += D[i, j]
    result = total / (len(ones) * len(twos))
    return(result)
def HierarchicalClustering(n, D):
    Clusters = [i for i in range(1, n + 1)]
    T = Graph()
    DCopy = np.copy(D)
    for i in range(n):
        T.newNode(i + 1)
    while len(Clusters) > 1:
        n = len(D)
        # find the two closest clusters
        flatD = [i for row in D for i in row]
        minDist = min([i for i in flatD if i != 0])
        index = flatD.index(minDist)
        i = index // n
        j = index % n
        clusterI = Clusters[i]
        clusterJ = Clusters[j]
        # add the new cluster to the graph
        new = str(clusterI) + " " + str(clusterJ)
        T.newEdge(new, clusterI)
        T.newEdge(new, clusterJ)
        # add row and col with updated distances
        entries = []
        for C in Clusters:
            ones = [int(x) - 1 for x in str(C).split(" ")]
            twos = [int(y) - 1 for y in new.split(" ")]
            dist = avgDist(DCopy, ones, twos)
            entries.append(dist)
        D = np.row_stack((D, entries))
        D = np.column_stack((D, entries + [0]))
        print(clusterI, clusterJ)
        # remove the old cluster rows and cols
        D = np.delete(D, [i, j], axis=0)
        D = np.delete(D, [i, j], axis=1)
        # update list of clusters
        Clusters.remove(clusterI)
        Clusters.remove(clusterJ)
        Clusters.append(new)
    # set the new root
    R = new
    T.setRoot(R)
    return(T)

with open('hw3/dataset_925937_7.txt') as f:
    n = int(f.readline().strip())
    lines = f.readlines()
    mat = []
    for line in lines:
        row = [float(i) for i in line.strip().split(' ')]
        mat.append(row)
dend = HierarchicalClustering(n, mat)