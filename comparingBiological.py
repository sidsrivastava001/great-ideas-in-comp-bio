import sys

sys.setrecursionlimit(10000)

"""
Algorithms for longest common subsequence
"""

"""
Backtrack through longest common subsequence DP table
Takes in two strings v and w
Returns backtrack matrix
"""
def LCSBackTrack(v, w):
    s = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    backtrack = [ [0]*(len(w)+1) for i in range(len(v)+1)]

    
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if(v[i-1] == w[j-1]):
                match = 1
            s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1]+match)
            if(s[i][j] == s[i-1][j]):
                backtrack[i][j] = "Down"
            elif(s[i][j] == s[i][j-1]):
                backtrack[i][j] = "Right"
            elif(s[i][j] == s[i-1][j-1]+match):
                backtrack[i][j] = "Diagonal"
    
    return backtrack
            
"""
Print out backtracked string through table
Takes in backtrack matrix, string v, and location (i, j) in matrix
Return concatenated longest common subsequence
"""
def OutputLCS(backtrack, v, i, j):
    if i == 0 or j == 0:
        return ""
    if(backtrack[i][j] == "Down"):
        return OutputLCS(backtrack, v, i - 1, j)
    if(backtrack[i][j] == "Right"):
        return OutputLCS(backtrack, v, i, j - 1)
    else:
        return OutputLCS(backtrack, v, i-1, j - 1)+v[i-1]

"""
Create DP table for longest common subsequence
Takes in 2 strings v, w, same bonus, mismatch penalty, and indel penalty
Returns backtrack DP matrix
"""
def LCSScores(v, w, same, mismatch, indel):
    s = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    backtrack = [ [0]*(len(w)+1) for i in range(len(v)+1)]
    for i in range(len(s)):
        s[i][0] = indel*i
    for i in range(len(s[0])):
        s[0][i] = indel*i
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if(v[i-1] == w[j-1]):
                
                match = same
            else:
                match = mismatch
            s[i][j] = max(s[i-1][j]+indel, s[i][j-1]+indel, s[i-1][j-1]+match)
            if(s[i][j] == s[i-1][j-1]+match):
                backtrack[i][j] = "Diagonal"
            elif(s[i][j] == s[i][j-1]+indel):
                backtrack[i][j] = "Right"
            else:
                backtrack[i][j] = "Down"
    print("SCORE", s[len(v)][len(w)])
    return backtrack
            
"""
Output LCS with scores
Takes in backtrack table, string v, location (i, j) in table
Return tuple of scores with concatenated LCS
"""
def outputWithScores(backtrack, v, i, j):
    if i == 0 or j == 0:
        print(i, j)
        return ["", ""]
    if(backtrack[i][j] == "Down"):
        x = outputWithScores(backtrack, v, i - 1, j)
        return [x[0]+"-", x[1]+v[i-1]]
    if(backtrack[i][j] == "Right"):
        x = outputWithScores(backtrack, v, i, j-1)
        return [x[0]+w[j-1], x[1]+"-"]
    else:
        x = outputWithScores(backtrack, v, i-1, j-1)
        return [x[0]+w[j-1], x[1]+v[i-1]]

if __name__ == '__main__':
    v = "GCGCATCGCTCATAGCTGTCCGGTCCATTGAAATGAACACTGTGACCCGGAATCCCGCGATGGAACCATCAACCCTCGATGTTTGTCTTTCGCCTGACTACGTTGGTATACCACCTTGACCCGATAGCCGGCGGCGTGTCCTAGACCAAATCATGCCCAACCTTTGGGACGTTCTAGAGCGATTAACATGGCCAAAGAGTTTTGTGCTAGGATAAAGAGCAAACTGAATTATCTGAGTCAGGGGCGTAAGCTTAGATGGGACACTATATGCAGAAAAAGATATACGGGGTGCTGAGCTCACACAGGAGGAACCCAGTCTTTGCTTGCAGACAACTTGGAAAAACGATGTAGACTGGTAGGAGTCAACGGATACGGGACACAAATAAGGAGCCGACTTGGTGAGCAATGGTCCGGCGGAACAGAATGAACGGGGTTGTTCAGCATTGGCATAGTTACTGGATCATCCGAGTCCTACGTGGACGAAGACGGGATTCGACTTGCTAAAATCCGTTATATCAAGCAGTTTATACCGATCTCGCTAACGGGACTTGTCCGGAGGCTGAGGGAATGAAACGTTTGCATGGGCAAGTCCCCACATGTACCTAAGTAACTGAGGAATCTAGCGGTGGCAGAATGATTTCTCCTACGCCGGTCACGTGTTCGAGTGCGGTATTAACTTTTGCATTACCAAAACTCAGATGCTACAAGGTGCTCTACCTAACTTCGGTTATACCATGATCTTCTCAACTGGGGTATGTGTTACCTGTAACTGTGGAACAAATGACTCTCCTGCGACGTGCAACAAACGGACATAGTCTAGGTGGTATTCGTGGCTTCTTGATGTCCCCACTCCTCAGTGAATTAGTCAATCGACAAGGACAATCTACCTCTACTCACTAC"
    w = "GCCTCTCAGTTTCGCTCCAGACAGCGTCCAGTCCATTGAAATGAAGTCGCCACACCCAATATTGTTAGACCCGTATATACGCTTTGTGGCTGATGGAACCATTAACCCTCGATGTTTGTCTTTCGCCTCTGAGATACGTTGGTATATCCACACGCGGGTTATGATAGCTGGTGGTGAGTCCTAGACCATGCCCACACGGCTGCCTTTGGGACGTTCTGATTAACATGGCGTCAAAGTGCAGCTAGGATCAAAGGCGTGCAGCCTGAGTCAGGGGAACTCGGTAAGACTGTTAGCTTAGAAGGGACATTATATGCAGAAAAATACACTGATGTATACGGAGTATTAACGCCTGAGCTCGGAGGAACCCAGTCTTTGCTTTCAGAAAAACGATGTAGCCTGGAAACCAACGGAAGAAGCGTACGGCACACAAATACCAAGGGAAGGTGCGTCTGACTAAGCAATGGTCCTAAATCGCAGTGTCTAGGGGTTGTTCAGCCTGCACGGTGGGCATAGTAAGATCATCCGAAGTCGGGATTCGACTGAAAATCCGTTATATTTATACCGATCTCGGGCCTCTAAAAACTATTGGGATGAAAAGAGAATGAGTCTTGACGATGCATGGGCACGTCCCTTTCTCGCTCCATATACCTAAGTAACAGAGGAATCTCGCGGTGGCAGAATCATGTTAAGTGGTCACGTGTTTCATATAATTTTTTTACCAAAACTCTGACTCGACGGTTCTCTCGCGCTCACTTCGGCGCTCGTTTTACTCCTCTCGCCCATGATCTTCTCAAATGTGTTAGTGCGCTGTGTCCCTGAACAAATGTCACACTCTCCTATCGCGCCGACGACACAGGTGGTGTGTAATATGCGTACGTGGTTTTCGTGGCTTCTCATCCTCCCGCTGAGTGAATTATTCAATCGACAAGGACATCGACGCATCTACCACTACTCACTGGC"
    backtrack = LCSScores(v, w, 1, -1, -5)
    x = outputWithScores(backtrack, v, len(v), len(w))
    f = open("output.txt", "w")
    f.write(x[1])
    f.write("\n")
    f.write(x[0])
    print(len(v), len(w))