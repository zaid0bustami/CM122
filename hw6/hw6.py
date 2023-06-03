#%% 5.4.3


#%% 5.6.6


#%% 5.6.10
import numpy as np
def ManhattanTourist(n, m, Down, Right):
    s = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        s[i][0] = s[i-1][0] + Down[i-1][0]
    for j in range(1, m + 1):
        s[0][j] = s[0][j-1] + Right[0][j-1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i-1][j] + Down[i-1][j], s[i][j-1] + Right[i][j-1])
    return(int(s[n][m]))


with open("hw6/dataset_925868_10.txt") as f:
    constants = [int(i) for i in f.readline().strip().split(" ")]
    n = constants[0]
    m = constants[1]
    Down = []
    while True:
        add = f.readline().strip()
        if add == "-":
            break
        Down.append(add.split(" "))
    Down = [[int(j) for j in i] for i in Down]
    Right = [[int(j) for j in i.strip().split(" ")] for i in f.readlines()]
    print(ManhattanTourist(n, m, Down, Right))

#%% 5.8.5
import numpy as np
import sys
sys.setrecursionlimit(10000)

def LCSBackTrack(v, w):
    s = np.empty((len(v) + 1, len(w) + 1), dtype=int)
    Backtrack = np.empty((len(v) + 1, len(w) + 1),  dtype=str)
    for i in range(len(v) + 1):
        s[i][0] = 0
    for j in range(len(w) + 1):
        s[0][j] = 0
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j]:
                Backtrack[i][j] = "|"
            elif s[i][j] == s[i][j-1]:
                Backtrack[i][j] = "-"
            elif s[i][j] == s[i-1][j-1] + match:
                Backtrack[i][j] = "\\"
    print(s)
    print(Backtrack)
    return(Backtrack)

def printBacktrack(Backtrack):
    for i in Backtrack:
        print(" ".join(i))
              
def OutputLCS(Backtrack, v, i, j):
    if i == 0 or j == 0:
        return("")
    if Backtrack[i][j] == "|":
        return(OutputLCS(Backtrack, v, i - 1, j))
    elif Backtrack[i][j] == "-":
        return(OutputLCS(Backtrack, v, i, j - 1))
    else:
        return(OutputLCS(Backtrack, v, i - 1, j - 1) + v[i - 1])
        print(i, j, Backtrack[i][j])


def LCS(v, w):
    Backtrack = LCSBackTrack(v, w)
    Output = OutputLCS(Backtrack, v, len(v), len(w))
    return(Output)

with open("hw6/dataset_925870_5.txt") as f:
    v = f.readline().strip()
    w = f.readline().strip()
    print(LCS(v, w))


#%% 5.10.3
import numpy as np

def AlignBacktrack(match, mismatch, indel, v, w):
    s = np.zeros((len(v) + 1, len(w) + 1), dtype=int)
    Backtrack = np.empty((len(v) + 1, len(w) + 1),  dtype=str)
    for i in range(1, len(v) + 1):
        s[i][0] = s[i-1][0] - indel
    for j in range(1, len(w) + 1):
        s[0][j] = s[0][j-1] - indel
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            diag = -mismatch
            if v[i-1] == w[j-1]:
                diag = match
            s[i][j] = max(s[i-1][j] - indel, s[i][j-1] - indel, s[i-1][j-1] + diag)
            if s[i][j] == s[i-1][j] - indel:
                Backtrack[i][j] = "|"
            elif s[i][j] == s[i][j-1] - indel:
                Backtrack[i][j] = "-"
            elif (s[i][j] == s[i-1][j-1] + match) or (s[i][j] == s[i-1][j-1] - mismatch):
                Backtrack[i][j] = "\\"
    return Backtrack, s[len(v)][len(w)]
                
def OutputAlignV(Backtrack, v, i, j):
    if i == 0 or j == 0:
        return("")
    if Backtrack[i][j] == "|":                          # deletion
        return(OutputAlignV(Backtrack, v, i - 1, j) + v[i-1])
    elif Backtrack[i][j] == "-":                        # insertion
        return(OutputAlignV(Backtrack, v, i, j - 1) + "-")
    else:                                               # (mis)match
        return(OutputAlignV(Backtrack, v, i - 1, j - 1) + v[i-1])

def OutputAlignW(Backtrack, w, i, j):
    if i == 0 or j == 0:
        return("")
    if Backtrack[i][j] == "|":                          # deletion
        return(OutputAlignW(Backtrack, w, i - 1, j) + "-")
    elif Backtrack[i][j] == "-":                        # insertion
        return(OutputAlignW(Backtrack, w, i, j - 1) + w[j-1])
    else:                                               # (mis)match
        return(OutputAlignW(Backtrack, w, i - 1, j - 1) + w[j-1])

def Align(match, mismatch, indel, v, w):
    Backtracks = AlignBacktrack(match, mismatch, indel, v, w)
    Backtrack = Backtracks[0]
    Score = Backtracks[1]
    OutputV = OutputAlignV(Backtrack, v, len(v), len(w))
    OutputW = OutputAlignW(Backtrack, w, len(v), len(w))
    return Score, OutputV, OutputW

with open("hw6/dataset_925872_3.txt") as f:
    scores = [int(i) for i in f.readline().strip().split(" ")]
    v = f.readline().strip()
    w = f.readline().strip()
    A = Align(scores[0], scores[1], scores[2], v, w)
    print(A[0])
    print(A[1])
    print(A[2])

#%% 5.11.3
def EditDistance(v, w):
    s = np.zeros((len(v) + 1, len(w) + 1), dtype=int)
    Backtrack = np.empty((len(v) + 1, len(w) + 1),  dtype=str)
    for i in range(1, len(v) + 1):
        s[i][0] = s[i-1][0] - 1
    for j in range(1, len(w) + 1):
        s[0][j] = s[0][j-1] - 1
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            diag = 1
            if v[i-1] == w[j-1]:
                diag = 0
            s[i][j] = max(s[i-1][j] - 1, s[i][j-1] - 1, s[i-1][j-1] - diag)
            if s[i][j] == s[i-1][j] - 1:
                Backtrack[i][j] = "|"
            elif s[i][j] == s[i][j-1] - 1:
                Backtrack[i][j] = "-"
            elif (s[i][j] == s[i-1][j-1]) or (s[i][j] == s[i-1][j-1] - 1):
                Backtrack[i][j] = "\\"
    return -s[len(v)][len(w)]    

with open("hw6/dataset_925873_3.txt") as f:
    v = f.readline().strip()
    w = f.readline().strip()
    print(EditDistance(v, w))

#%% 5.13.12 Extra Credit
