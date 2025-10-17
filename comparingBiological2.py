import sys
import math
sys.setrecursionlimit(10000)

"""
More longest common subsequence algorithms with different penalties
Also includes HMM implmentation for multiple protein sequence alignment
"""

"""
Create DP table for longest common subsequence
Takes in 2 strings v, w, mismatch penalty, and indel penalty
Returns backtrack table from scores
"""
def LCSScores(v, w, mismatch, indel):
    ma = 0
    maxpair = [0, 0]
    s = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    backtrack = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    for i in range(len(s)):
        s[i][0] = indel*i
    for i in range(len(s[0])):
        s[0][i] = indel*i
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = mismatch[(v[i-1], w[j-1])]
            s[i][j] = max(s[i-1][j]+indel, s[i][j-1]+indel, s[i-1][j-1]+match, 0)
            if(s[i][j] > ma):
                ma = s[i][j]
                maxpair = [i, j]
            if(s[i][j] == s[i-1][j-1]+match):
                backtrack[i][j] = "Diagonal"
            elif(s[i][j] == s[i][j-1]+indel):
                backtrack[i][j] = "Right"
            elif(s[i][j] == s[i-1][j]+indel):
                backtrack[i][j] = "Down"
            else:
                backtrack[i][j] = "Back"
    print("SCORE", s[maxpair[0]][maxpair[1]])
    return backtrack, maxpair

"""
Print out backtracked string through table with scores
Takes in backtrack table, string v, and position i, j
Returns backtracked scores from table
"""
def outputWithScores(backtrack, v, i, j):
    if i == 0 or j == 0:
        print(i, j)
        return ["", ""]
    if(backtrack[i][j] == "Down"):
        x = outputWithScores(backtrack, v, i - 1, j)
        return [x[0]+"-", x[1]+v[i-1]]
    elif(backtrack[i][j] == "Right"):
        x = outputWithScores(backtrack, v, i, j-1)
        return [x[0]+w[j-1], x[1]+"-"]
    elif(backtrack[i][j] == "Diagonal"):
        x = outputWithScores(backtrack, v, i-1, j-1)
        return [x[0]+w[j-1], x[1]+v[i-1]]
    else:
        return outputWithScores(backtrack, v, 0, 0)



"""
Create DP table for longest common subsequence with gap opening and extension penalties
Takes in strings v and w, mismatch penalty, same nucleotide bonus, gap_open penalty, and gap_extension penalty
Returns backtrack table with scores
"""
def LCSScoresAffine(v, w, same, mismatch, gap_open, gap_extension):
    s = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    upper = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    lower = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    backtrack = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    for i in range(len(v)+1):
        for j in range(len(w)+1):
            if(i == 0 and j != 0):
                s[i][j] = gap_open+(j-1)*gap_extension
                lower[i][j] = -999999
            if(i != 0 and j == 0):
                s[i][j] = gap_open+(i-1)*gap_extension
                upper[i][j] = -999999
            if(i ==0 and j == 0):
                upper[i][j] = -999999
                lower[i][j] = -999999

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if(v[i-1] == w[j-1]):
                match = same
            else:
                match = mismatch
            lower[i][j] = max(lower[i-1][j]+gap_extension, s[i-1][j]+gap_open)
            upper[i][j] = max(upper[i][j-1]+gap_extension, s[i][j-1]+gap_open)
            s[i][j] = max(lower[i][j], upper[i][j], s[i-1][j-1]+match)
            if(s[i][j] == s[i-1][j-1]+match):
                backtrack[i][j] = "Diagonal"
            elif(s[i][j] == upper[i][j]):
                backtrack[i][j] = "Right"
            else:
                backtrack[i][j] = "Down"
    for i in s:
        print(i)
    for i in upper:
        print(i)
    for i in lower:
        print(i)
    print("SCORE", max(s[len(v)][len(w)], upper[len(v)][len(w)], lower[len(v)][len(w)]))
    print("Other score", s[len(v)][len(w)])
    return backtrack
            
"""
Output LCS scores with gap opening/extension
Takes in backtrack table, string v, and current position i, j
Returns backtracked scores
"""
def outputWithScoresAffine(backtrack, v, i, j):
    if i == 0 and j == 0:
        print(i, j)
        return ["", ""]
    if(backtrack[i][j] == "Down"):
        x = outputWithScoresAffine(backtrack, v, i - 1, j)
        return [x[0]+"-", x[1]+v[i-1]]
    if(backtrack[i][j] == "Right"):
        x = outputWithScoresAffine(backtrack, v, i, j-1)
        return [x[0]+w[j-1], x[1]+"-"]
    else:
        x = outputWithScoresAffine(backtrack, v, i-1, j-1)
        return [x[0]+w[j-1], x[1]+v[i-1]]


"""
Create 3D DP Table for multiple string longest subsequence
Takes in 3 strings v, w, and x
Returns backtrack table with scores
"""
def LCSBackTrackMultiple(v, w, x):
    s = [[[0 for k in range(len(x)+1)] for j in range(len(w)+1)] for i in range(len(v)+1)]
    backtrack = [[[7 for k in range(len(x)+1)] for j in range(len(w)+1)] for i in range(len(v)+1)]

    for i in range(len(v)+1):
        for j in range(len(w)+1):
            for k in range(len(x)+1):
                if(i == 0 and j == 0):
                    backtrack[i][j][k] = 3
                elif(i == 0 and k == 0):
                    backtrack[i][j][k] = 2
                elif(j == 0 and k == 0):
                    backtrack[i][j][k] = 1
                elif(i == 0):
                    backtrack[i][j][k] = 6
                elif(j == 0):
                    backtrack[i][j][k] = 5
                elif(k == 0):
                    backtrack[i][j][k] = 4

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            for k in range(1, len(x)+1):
                match = 0
                if(v[i-1] == w[j-1] and w[j-1] == x[k-1]):
                    match = 1
                s[i][j][k] = max(s[i-1][j][k], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k], s[i][j-1][k-1], s[i][j][k-1], s[i-1][j-1][k-1]+match)
                if(s[i][j][k] == s[i-1][j][k]):
                    backtrack[i][j][k] = 1
                elif(s[i][j][k] == s[i][j-1][k]):
                    backtrack[i][j][k] = 2
                elif(s[i][j][k] == s[i][j][k-1]):
                    backtrack[i][j][k] = 3
                elif(s[i][j][k] == s[i-1][j-1][k]):
                    backtrack[i][j][k] = 4
                elif(s[i][j][k] == s[i-1][j][k-1]):
                    backtrack[i][j][k] = 5
                elif(s[i][j][k] == s[i][j-1][k-1]):
                    backtrack[i][j][k] = 6
                elif(s[i][j][k] == s[i-1][j-1][k-1]+match):
                    backtrack[i][j][k] = 7
    print("SCORE: ", s[len(v)][len(w)][len(x)])
    for i in backtrack:
        print(i)
    return backtrack


"""
Output LCS for multiple strings with 3D DP
Takes in backtrack matrix, strings v, w, and x, and position in 3D DP table (i, j, k)
Returns concatenated string of best alignment
"""         
def OutputLCSMultiple(backtrack, v, w, x, i, j, k):
    if i == 0 and j == 0 and k == 0:
        return ["", "", ""]
    if(backtrack[i][j][k] == 1):
        y = OutputLCSMultiple(backtrack, v, w, x, i - 1, j, k)
        print(y[0]+v[i-1], y[1]+"-", y[2]+"-")
        return [y[0]+v[i-1], y[1]+"-", y[2]+"-"]
    elif(backtrack[i][j][k] == 2):
        y = OutputLCSMultiple(backtrack, v, w, x, i, j-1, k)
        print(y[0]+"-", y[1]+w[j-1], y[2]+"-")
        return [y[0]+"-", y[1]+w[j-1], y[2]+"-"]
    elif(backtrack[i][j][k] == 3):
        y = OutputLCSMultiple(backtrack, v, w, x, i, j, k-1)
        print(y[0]+"-", y[1]+"-", y[2]+x[k-1])
        return [y[0]+"-", y[1]+"-", y[2]+x[k-1]]
    elif(backtrack[i][j][k] == 4):
        y = OutputLCSMultiple(backtrack, v, w, x, i-1, j-1, k)
        print(y[0]+v[i-1], y[1]+w[j-1], y[2]+"-")
        return [y[0]+v[i-1], y[1]+w[j-1], y[2]+"-"]
    elif(backtrack[i][j][k] == 5):
        y = OutputLCSMultiple(backtrack, v, w, x, i-1, j, k-1)
        print(y[0]+v[i-1], y[1]+"-", y[2]+x[k-1])
        return [y[0]+v[i-1], y[1]+"-", y[2]+x[k-1]]
    elif(backtrack[i][j][k] == 6):
        y = OutputLCSMultiple(backtrack, v, w, x, i, j-1, k-1)
        print(y[0]+"-", y[1]+w[j-1], y[2]+x[k-1])
        return [y[0]+"-", y[1]+w[j-1], y[2]+x[k-1]]
    elif(backtrack[i][j][k] == 7):
        y = OutputLCSMultiple(backtrack, v, w, x, i-1, j-1, k-1)
        print(y[0]+v[i-1], y[1]+w[j-1], y[2]+x[k-1])
        return [y[0]+v[i-1], y[1]+w[j-1], y[2]+x[k-1]]

"""
Hidden Markov Model for multiple protein sequence alignment
Takes in transmission probabilities, emmision probabilities, list of protein strings, and numer of states
Returns most probable aligned protein sequence
"""
def HMM(transmission, emission, string, states):
    weights = [ [0]*(len(states)) for i in range(len(string))]
    backtrack = [ [0]*(len(states)) for i in range(len(string))]

    for i in range(len(states)):
        weights[0][i] = 1/len(states)
    for i in weights:
        print(i)
    for i in range(1, len((string))):
        print(i)
        for j in range(len(states)):
            max = -100000000
            currState = states[j]
            for k in range(len(states)):
                prevState = states[k]
                if(math.log(transmission[prevState][currState]*emission[currState][string[i]])+(weights[i-1][k]) > max):
                    max = math.log(transmission[prevState][currState]*emission[currState][string[i]])+(weights[i-1][k]) 
                    weights[i][j] = max
                    backtrack[i][j] = k
                    print(max, i, j, k)
    max = -10000
    maxindex = 0
    for i in range(len(states)):
        if(weights[len(string)-1][i] > max):
            max = weights[len(string)-1][i]
            maxindex = i
    ret = states[maxindex]
    prev = maxindex
    for i in backtrack:
        print(i)
    for i in weights:
        print(i)
    for i in range(len(string)-1, 0, -1):
        print(states[backtrack[i][prev]])
        ret = states[backtrack[i][prev]]+ret
        prev = backtrack[i][prev]
    return ret

