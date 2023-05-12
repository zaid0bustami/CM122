#%% Libraries
from sklearn.linear_model import LogisticRegression

# Functions
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
    return maxProb[1], maxProb[0]

def kmerComp(Txt, k):
    matches = dict()
    for kmer in kmers(k):
        matches[kmer] = 0
        for i in range(len(Txt) - k + 1):
            t = Txt[i:i+k]
            if t == kmer:
                matches[kmer] += 1
    return(matches)
            
def kmers(k):
    bases = {"A", "T", "C", "G"}
    if k <= 1:
        return(bases)
    else:
        result = []
        for b in bases:
            for rest in kmers(k - 1):
                result.append(b + rest)
        return(result)

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

 #%% Calculating Predictions
bound = readsFASTA("data/PROJECT2B/bound.fasta")
notbound = readsFASTA("data/PROJECT2B/notbound.fasta")
test = readsFASTA("data/PROJECT2B/test.fasta")
k = 3

print("Pre-processing Training Data")
train = [[i[1], 1] for i in bound] + [[i[1], 0] for i in notbound]
X_train = [[count for kmer, count in kmerComp(i[0], k).items()] for i in train]
y_train = [i[1] for i in train]

print("Pre-processing Test Data")
X_test = [[count for kmer, count in kmerComp(i[1], k).items()] for i in test]

print("Running Logistic Regression")
log_reg = LogisticRegression(solver='lbfgs', multi_class='auto', max_iter=5000)
log_reg.fit(X_train, y_train)
y_pred = log_reg.predict_proba(X_test)[:,1]

print("Returning Predictions")
predictions = [(test[i][0], y_pred[i]) for i in range(len(test))]

#%% Filtering for Top 2000
sortedPredictions = sorted(predictions, key = lambda x: x[1], reverse=True)

finalPredictions = [i[0][1:] for i in sortedPredictions[0:2000]]

with open('predictions.csv', 'w') as o:
    for i in finalPredictions:
        o.write(i + "\n")
