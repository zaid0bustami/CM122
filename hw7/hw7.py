#%% 10.5.8
def probPath(path, states, transition):
    index = dict()
    result = 1 / len(states)
    for i in range(len(states)):
        add = {states[i] : i}
        index.update(add)
    for i in range(len(path) - 1):
        state1 = path[i]
        state2 = path[i+1]
        prob = transition[index[state1]][index[state2]]
        print(state1, state2, prob)
        result *= prob
    return result

with open("hw7/dataset_925971_8.txt") as f:
    path = f.readline().strip()
    f.readline()
    states = f.readline().strip().split(" ")
    f.readline()
    transition = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]
    prob = probPath(path, states, transition)
    print(prob)
        

#%% 10.5.10
def probEmissionsGivenPath(string, alphabet, path, states, emission):
    alphaIndex = dict()
    stateIndex = dict()
    result = 1
    for i in range(len(alphabet)):
        add = {alphabet[i] : i}
        alphaIndex.update(add)
    for i in range(len(states)):
        add = {states[i] : i}
        stateIndex.update(add)
    for i in range(len(string)):
        emit = string[i]
        state = path[i]
        prob = emission[stateIndex[state]][alphaIndex[emit]]
        print(emit, state, prob)
        result *= prob
    return(result)


with open("hw7/dataset_925971_10.txt") as f:
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split(" ")
    f.readline()
    path = f.readline().strip()
    f.readline()
    states = f.readline().strip().split(" ")
    f.readline()
    emission = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]
    prob = probEmissionsGivenPath(string, alphabet, path, states, emission)
    print(prob)

#%% 10.6.7
import numpy as np

def Weight(emit, statePrev, stateCur, 
           alphaIndex, stateIndex, transition, emission, 
           probInitialTransition = None):
    if statePrev != None:
        probTransition = transition[stateIndex[statePrev]][stateIndex[stateCur]]
    else:
        probTransition = probInitialTransition
    probEmission = emission[stateIndex[stateCur]][alphaIndex[emit]]
    result = probTransition * probEmission
    return result
    
def ViterbiDP(string, states, alphaIndex, stateIndex, transition, emission):
    # initializing
    s = np.zeros((len(states), len(string)), dtype=float)
    backtrack = np.zeros((len(states), len(string) + 1), dtype=str)
    for i in range(len(states)):
        initialWeight = Weight(string[0], None, states[i], 
                               alphaIndex, stateIndex, transition, emission,
                               probInitialTransition = 1 / len(states))
        s[i][0] = initialWeight
    # dynamic programming and traversal
    for j in range(1, len(string)):
        for i in range(len(states)):
            emit = string[j]
            stateCur = states[i]
            weights = [Weight(emit, statePrev, stateCur, 
                             alphaIndex, stateIndex, transition, emission)
                       for statePrev in states]
            allWeights = [weights[k] * s[k][j-1]
                          for k in range(len(states))]
            s[i][j] = max(allWeights)
            statePrev = states[allWeights.index(s[i][j])]
            backtrack[i][j] = statePrev
    # adding final entry to backtrack
    lastWeights = [i[-1] for i in s]
    lastStateIndex = lastWeights.index(max(lastWeights))
    lastState = states[lastStateIndex]
    for i in range(len(states)):
        backtrack[i][-1] = lastState
    return backtrack

def ViterbiPath(backtrack, states, i, j):
    if j == 0:
        return("")
    prevState = backtrack[i][j]
    prevIndex = states.index(prevState)
    return(ViterbiPath(backtrack, states, prevIndex, j - 1) + prevState)

def Viterbi(string, alphabet, states, transition, emission):
    alphaIndex = dict()
    stateIndex = dict()
    for i in range(len(alphabet)):
        add = {alphabet[i] : i}
        alphaIndex.update(add)
    for i in range(len(states)):
        add = {states[i] : i}
        stateIndex.update(add)
    backtrack = ViterbiDP(string, states, 
                          alphaIndex, stateIndex, transition, emission)
    path = ViterbiPath(backtrack, states, 0, len(string))
    return(path)

with open("hw7/dataset_925972_7.txt") as f:
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split(" ")
    f.readline()
    states = f.readline().strip().split(" ")
    f.readline()
    f.readline()
    transition = []
    for i in range(4):
        transition.append([float(j) for j in f.readline().strip().split("\t")[1:]])
    f.readline()
    emission = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]
    path = Viterbi(string, alphabet, states, transition, emission)
    print(path)

#%% 10.7.4
def ForwardDP(string, states, alphaIndex, stateIndex, transition, emission):
    # initializing
    s = np.zeros((len(states), len(string)), dtype=float)
    for i in range(len(states)):
        initialWeight = Weight(string[0], None, states[i], 
                               alphaIndex, stateIndex, transition, emission,
                               probInitialTransition = 1 / len(states))
        s[i][0] = initialWeight
    # dynamic programming and traversal
    for j in range(1, len(string)):
        for i in range(len(states)):
            emit = string[j]
            stateCur = states[i]
            weights = [Weight(emit, statePrev, stateCur, 
                             alphaIndex, stateIndex, transition, emission)
                       for statePrev in states]
            allWeights = [weights[k] * s[k][j-1]
                          for k in range(len(states))]
            s[i][j] = sum(allWeights)
    # returning the forward result
    lastWeights = [i[-1] for i in s]
    return sum(lastWeights)

def Forward(string, alphabet, states, transition, emission):
    alphaIndex = dict()
    stateIndex = dict()
    for i in range(len(alphabet)):
        add = {alphabet[i] : i}
        alphaIndex.update(add)
    for i in range(len(states)):
        add = {states[i] : i}
        stateIndex.update(add)
    result = ForwardDP(string, states, 
                       alphaIndex, stateIndex, transition, emission)
    return result

with open("hw7/dataset_925973_4.txt") as f:
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split(" ")
    f.readline()
    states = f.readline().strip().split(" ")
    f.readline()
    f.readline()
    transition = []
    for i in range(4):
        transition.append([float(j) for j in f.readline().strip().split("\t")[1:]])
    f.readline()
    emission = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]
    result = Forward(string, alphabet, states, transition, emission)
    print(result)

#%% 10.11.4
def SupervisedHMMLearning(string, alphabet, path, states):
    # initialize matrices
    alphaIndex = dict()
    stateIndex = dict()
    for i in range(len(alphabet)):
        add = {alphabet[i] : i}
        alphaIndex.update(add)
    for i in range(len(states)):
        add = {states[i] : i}
        stateIndex.update(add)
    transition = np.zeros((len(states), len(states)), dtype=float)
    emission = np.zeros((len(states), len(alphabet)), dtype=float)
    # fill in emission matrix
    nStates = []
    for state in states:
        filteredString = []
        for i in range(len(path)):
            s = path[i]
            x = string[i]
            if s == state:
                filteredString.append(x)
        pseudo = 0
        if len(filteredString) == 0:
            pseudo = 1
        nState = [s for s in path].count(state) + pseudo * len(alphabet)
        nStates.append(nState)
        for emit in alphabet:
            nEmit = filteredString.count(emit) + pseudo
            emission[stateIndex[state]][alphaIndex[emit]] = nEmit / nState      
    # fill in transition matrix
    for i in range(len(path) - 1):
        state1 = path[i]
        state2 = path[i + 1]
        transition[stateIndex[state1]][stateIndex[state2]] += 1
    for i in range(len(transition)):
        row = transition[i]
        if sum(row) == 0:
            row = [r + 1 for r in row]
        denom = sum(row)
        transition[i] = [r / denom for r in row]
    return transition, emission

def Display(result):
    transition = result[0]
    emission = result[1]
    for s in states:
        print("\t", end = s)
    print()
    for i in range(len(states)):
        s = states[i]
        print(s, end = "")
        for j in transition[i]:
            print("\t", end = str(round(j, 3)))
        print()
    print("--------")
    for x in alphabet:
        print("\t", end = x)
    print()
    for i in range(len(states)):
        s = states[i]
        print(s, end = "")
        for j in emission[i]:
            print("\t", end = str(round(j, 3)))
        print()

with open("hw7/dataset_925977_4.txt") as f:
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split(" ")
    f.readline()
    path = f.readline().strip()
    f.readline()
    states = f.readline().strip().split(" ")
    result = SupervisedHMMLearning(string, alphabet, path, states)
    Display(result)

#%% 10.12.5
def ForwardMat(string, states, alphaIndex, stateIndex, transition, emission):
    # initializing
    s = np.zeros((len(states), len(string)), dtype=float)
    for i in range(len(states)):
        initialWeight = Weight(string[0], None, states[i], 
                               alphaIndex, stateIndex, transition, emission,
                               probInitialTransition = 1 / len(states))
        s[i][0] = initialWeight
    # dynamic programming and traversal
    for j in range(1, len(string)):
        for i in range(len(states)):
            emit = string[j]
            stateCur = states[i]
            weights = [Weight(emit, statePrev, stateCur, 
                             alphaIndex, stateIndex, transition, emission)
                       for statePrev in states]
            allWeights = [weights[k] * s[k][j-1]
                          for k in range(len(states))]
            s[i][j] = sum(allWeights)
    return s

def BackwardMat(string, states, alphaIndex, stateIndex, transition, emission):
    # initializing
    s = np.zeros((len(states), len(string)), dtype=float)
    for i in range(len(states)):
        initialWeight = 1
        s[i][-1] = initialWeight
    # dynamic programming and traversal
    for j in reversed(range(len(string) - 1)):
        for i in range(len(states)):
            emit = string[j + 1]
            stateCur = states[i]
            weights = [Weight(emit, stateCur, statePrev, 
                             alphaIndex, stateIndex, transition, emission)
                       for statePrev in states]
            allWeights = [weights[k] * s[k][j+1]
                          for k in range(len(states))]
            s[i][j] = sum(allWeights)
    return s

def ForwardBackward(string, alphabet, states, transition, emission):
    alphaIndex = dict()
    stateIndex = dict()
    for j in range(len(alphabet)):
        add = {alphabet[j] : j}
        alphaIndex.update(add)
    for j in range(len(states)):
        add = {states[j] : j}
        stateIndex.update(add)
    result = np.zeros((len(states), len(string)), dtype=float)
    forwardMat = ForwardMat(string, states, alphaIndex, stateIndex, transition, emission)
    backwardMat = BackwardMat(string, states, alphaIndex, stateIndex, transition, emission)
    forwardSink = sum([i[-1] for i in forwardMat])
    result = np.multiply(forwardMat, backwardMat) / forwardSink       
    return(result)         

with open("hw7/dataset_925978_5.txt") as f:
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split(" ")
    f.readline()
    states = f.readline().strip().split(" ")
    f.readline()
    f.readline()
    transition = []
    for i in range(len(states)):
        transition.append([float(j) for j in f.readline().strip().split("\t")[1:]])
    f.readline()
    emission = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]

    SoftDecoding = ForwardBackward(string, alphabet, states, transition, emission)
    result = np.transpose(SoftDecoding)
    print("\t".join(states))
    test = ["\t".join([str(round(c, 3)) for c in r]) for r in result]
    for r in test:
        print(r)

#%% 10.13.5
def BaumWelch(j, string, alphabet, states, transition, emission):
    # initializing
    alphaIndex = dict()
    stateIndex = dict()
    for i in range(len(alphabet)):
        add = {alphabet[i] : i}
        alphaIndex.update(add)
    for i in range(len(states)):
        add = {states[i] : i}
        stateIndex.update(add)
    M = [transition, emission]
    # iterating
    for iteration in range(j):
        # calculating forward and backward variables
        forwardMat = ForwardMat(string, states, alphaIndex, stateIndex, M[0], M[1])
        backwardMat = BackwardMat(string, states, alphaIndex, stateIndex, M[0], M[1])
        forwardSink = sum([i[-1] for i in forwardMat])
        # estimating parameters
        PiStar = np.multiply(forwardMat, backwardMat) / forwardSink
        PiStarStar = np.zeros((len(states)**2, len(string) - 1), dtype=float)
        for l in range(len(string) - 1):
            r = 0
            for state1 in states:
                for state2 in states:
                    forward = forwardMat[stateIndex[state1]][l]
                    weight = Weight(string[l], state1, state2, 
                                    alphaIndex, stateIndex, M[0], M[1])
                    backward = backwardMat[stateIndex[state2]][l + 1]
                    PiStarStar[r][l] = (forward * weight * backward) / forwardSink
                    r += 1
        # updating transition
        r = 0
        for state1 in states:
            for state2 in states:
                newProb = sum(PiStarStar[r])
                M[0][stateIndex[state1]][stateIndex[state2]] = newProb
                r += 1
        M[0] = M[0] / np.sum(M[0])
        # updating emission
        sums = dict()
        for s in states:
            innersum = dict()
            for x in alphabet:
                innersum.update({x : 0})
            sums.update({s : innersum})
        for k in range(len(states)):
            for m in range(len(string)):
                emit = string[m]
                sums[states[k]][emit] += PiStar[k][m]
        for s in states:
            for x in alphabet:
                M[1][stateIndex[s]][alphaIndex[x]] = sums[s][x]
        M[1] = M[1] / np.sum(M[1])
    return M

with open("hw7/dataset_925979_5.txt") as f:
    j = int(f.readline().strip())
    f.readline()
    string = f.readline().strip()
    f.readline()
    alphabet = f.readline().strip().split("\t")
    f.readline()
    states = f.readline().strip().split("\t")
    f.readline()
    f.readline()
    transition = []
    for i in range(len(states)):
        transition.append([float(j) for j in f.readline().strip().split("\t")[1:]])
    f.readline()
    emission = [[float(j) for j in i.strip().split("\t")[1:]] for i in f.readlines()[1:]]
    
    j = 10
    string = "xzyyzyzyxy"
    alphabet = ["x",	"y",	"z"]
    states = ["A",	"B"]
    transition = [[0.019,	0.981], 
                  [0.668,	0.332]]
    emission = [[0.175,	0.003,	0.821],
                [0.196,	0.512,	0.293]]
    
    result = BaumWelch(j, string, alphabet, states, transition, emission)
    Display(result)
    
#%%
from hmmlearn import hmm
import numpy as np

# Initialize the HMM model
model = hmm.MultinomialHMM(n_components=3, n_iter=100)

# Training data (observed sequence)
X = np.array([[0, 1, 0, 2]])

# Fit the model to the training data using the Baum-Welch algorithm
model.fit(X)

# Get the estimated model parameters
print("Estimated initial state probabilities:", model.startprob_)
print("Estimated transition matrix:", model.transmat_)
print("Estimated emission matrix:", model.emissionprob_)
