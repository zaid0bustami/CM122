#%% Libraries
from Bio import SeqIO as io
import pandas as pd

#%% Classes and Functions
class BWT:
    def __init__(self, Text):
        #compute suffix array and BWT
        print("\tCycles Array")
        Array = []
        for i in range(len(Text)):
            # entry = (i, Text[i:len(Text)] + Text[0:i])
            entry = (i, Text[i] + Text[i - 1])
            Array.append(entry)
        print("\tSorting")
        sortedArray = sorted(Array, key = lambda x: x[1])
        bwtList = [i[1][1] for i in sortedArray]
        self.SuffixArray = [i[0] for i in sortedArray]
        self.BWT = "".join(bwtList)
                
        # compute lasttofirst array
        print("\tLast-to-First Array")
        rightArray = [bwtList[i] + str(i) for i in range(len(bwtList))]
        leftArray = sorted(rightArray, key = lambda x: x[0])
        self.LF = [int(i[1:]) for i in leftArray]
        
    def Hamming(self, t, u):
        h = 0
        for i in range(len(t)):
            if t[i] != u[i]:
                h += 1
        return(h)

    def Neighbors(self, t, d):
        bases = {"A", "C", "G", "T"}
        if d == 0:
            return({t})
        if len(t) == 1:
            return(bases)
        neighborhood = set()
        suffixT = t[1:len(t)]
        suffixNeighbors = self.Neighbors(suffixT, d)
        for s in suffixNeighbors:
            if self.Hamming(suffixT, s) < d:
                for b in bases:
                    neighborhood.add(b + s)
            else:
                neighborhood.add(t[0] + s)
        return(neighborhood)
    
    def match(self, Pattern, d):
        result = []
        neighborhood = self.Neighbors(Pattern, d)
        result = [(n, self.matchHelper(n)) for n in neighborhood]
        print("matched")
        # for n in neighborhood:
        #     matches = self.matchHelper(n)
        #     result.append((n, matches))
        return(result)

    def matchHelper(self, Pattern):
        top = 0
        bottom = len(self.BWT) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                if len(Pattern) % 100 == 0:
                    print(".", end="")
                symbol = Pattern[-1]
                Pattern = Pattern[0:-1]
                topToBottom = self.BWT[top : bottom + 1]
                if symbol in topToBottom:
                    topIndex = top + topToBottom.index(symbol)
                    for i in range(len(topToBottom) - 1, -1, -1):
                        if topToBottom[i] == symbol:
                            bottomIndex = top + i
                            break
                    top = self.LF[topIndex]
                    bottom = self.LF[bottomIndex]
                else:
                    return -1
            else:
                return(self.SuffixArray[top:bottom + 1])

def FASTA(path):
    ID = {}
    for i in io.parse(path, "fasta"):
        ID[i.id] = i.seq
    return(ID)
def genomeFASTA(path):
    result = str(FASTA(path)['genome'])
    return(result)
def readsFASTA(path):
    with open(path) as f:
        reads = [read.strip() for read in f.readlines() if read[0] != '>']
        return(reads)
    
#%% Inputs
# sampleGenome = genomeFASTA('Data\\sample_1000\\sample_1000_reference_genome.fasta')
# sampleReads = readsFASTA('Data\\sample_1000\\sample_1000_no_error_single_reads.fasta')
# sampleReadsError = readsFASTA('Data\\sample_1000\\sample_1000_with_error_single_reads.fasta')
# sampleReadsPaired = readsFASTA('Data\\sample_1000\\sample_1000_no_error_paired_reads.fasta')
# sampleReadsErrorPaired = readsFASTA('Data\\sample_1000\\sample_1000_with_error_paired_reads.fasta')

genome = genomeFASTA('Data\\project1a_10000_reference_genome.fasta')
reads = readsFASTA('Data\\project1a_10000_with_error_paired_reads.fasta')

genome = genomeFASTA('Data\\project1b_1000000_reference_genome.fasta')
reads = readsFASTA('Data\\project1b_1000000_with_error_paired_reads.fasta')

#%% Read Matching
print('Storing BWT')
GENOME = genome
READS = reads
genomeBWT = BWT(GENOME + "$")

# create a list of dictionaries, where each contains keys as read thirds and 
# values as their indices
print('Mapping Reads')
readMap = []
for read in READS:
    n = len(read)
    readThirds = [read[:n//3], read[n//3:2*n//3], read[2*n//3:]]
    if(read != "".join(readThirds)):
        raise Exception
    for readThird in readThirds:
        matches = genomeBWT.match(readThird, 0)
        readMap += matches
    #%%    
print('Building Dataframe')
df = pd.DataFrame(columns=['base', 'reference', 'position'])
for i in range(len(readMap)):
    # print(".", end = "")
    match = readMap[i]
    read = match[0]
    start = match[1]
    length = len(read)
    if start != -1:                     # read counts for complete matches
        bases = list(read)
        positions = [i for i in range(start[0], start[0] + length)]
        references = [GENOME[i] for i in positions]
        ef = pd.DataFrame({'base': bases, 
                            'reference' : references,
                            'position': positions})
        df = pd.concat([df, ef])
    else:                               # read counts for mismatches
        whichThird = i % 3 
        matches = readMap[ (i - whichThird) : (i - whichThird + 3) ]
        reads = [i[0] for i in matches]
                
        nextThird = (whichThird + 1) % 3
        nextMatch = matches[nextThird]
        nextStart = nextMatch[1]
        if nextStart == -1:
            nextThird = (nextThird + 1) % 3
            nextMatch = matches[nextThird]
            nextStart = nextMatch[1]
        if nextStart == -1:
            continue
        
        nextRead = nextMatch[0]
        nextStart = nextStart[0] # assuming there's only one match per read
        shift = nextThird - whichThird
        if shift == 1:
            actualStart = nextStart - length
        elif shift == 2:
            actualStart = nextStart - length - len(reads[nextThird - 1])
        elif shift == -1:
            actualStart = nextStart + len(nextRead)
        elif shift == -2:
            actualStart = nextStart + len(nextRead) + len(reads[nextThird + 1])
        
        # print(whichThird, nextThird, shift, nextStart, actualStart)
        
        bases = list(read)
        positions = [i for i in range(actualStart, actualStart + length)]
        references = [GENOME[i] for i in positions]
        
        ef = pd.DataFrame({'base': bases, 
                            'reference' : references,
                            'position': positions})
        df = pd.concat([df, ef])
        
df = df.sort_values('position')
# df = (df.groupby(['base', 'reference', 'position', 'category'])
#       .size()
#       .rename('count')
#       .reset_index(drop=False)
#       .sort_values('position')
#       )

#%% Find Substitutions
# mismatchList = []
with open('predictions.txt', 'w') as f:
    positions = df.position.unique()
    prevPos = -5
    for i in range(len(positions)):
        pos = positions[i]
        dfSlice = df[df.position == pos]
        propMismatch = sum(dfSlice.base != dfSlice.reference) / len(dfSlice)
        if (propMismatch > 0.8): # we have a mutation
            if abs(pos - prevPos) > 2: # we have a substitution
                output = ">S" + str(dfSlice.iat[0, 2]) + " " + dfSlice.iat[0, 1] + " " + dfSlice.iat[0, 0]
                print(".", end="")
                print(output, file=f)
            # mismatchList.append(dfSlice)
            prevPos = pos

# mv ./predictions.txt ./predictions.csv