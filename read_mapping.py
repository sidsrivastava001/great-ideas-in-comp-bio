import math

"""
Functions to map reads using burrows-wheel transform, suffix arrays/trees/tries
"""

"""
Create a burrows wheeler transform matrix of a string
Takes in string
Returns BWT of string with most repeatable elements
"""
def BWT(str):
    matrix = []
    for i in range(len(str)):
        matrix.append(str[i:len(str)]+str[0:i])
    matrix.sort()
    bwt = ""
    for i in matrix:
        bwt+=(i[len(str)-1])
    return bwt

"""
Reverse burrows wheeler transform given a matrix
Takes in bwt matrix
Return original string
"""
def ReverseBWT(bwt):
    first = ''.join(sorted(bwt))
    matrix = []
    for i in first:
        matrix.append(i)
    for i in range(len(bwt)-1):
        for j in range(len(bwt)):
            matrix[j] = bwt[j]+matrix[j]
            matrix.sort()
    for i in range(1, len(bwt)):
        print(matrix[0][i], end = "")
    print(matrix[0][0])


"""
Create a suffix array from a string
Takes in string
Returns sorted dict of all suffixes
"""
def suffixArray(str):
    d = {}
    for i in range(len(str)):
        d[str[i:len(str)]] = i
    myKeys = list(d.keys())
    myKeys.sort()
    sorted_dict = {i: d[i] for i in myKeys}
    return sorted_dict


"""
Create a suffix trie from a string
Takes in string
Return list of nodes, edges dict, and 
"""
def suffixTrie(str):
    edges = {}
    nodes = [0]
    edge_thing = {}
    for i in range(len(str)):
        current = nodes[0]
        for j in range(i, len(str)):
            flag = False
            if(current in edge_thing):
                for k in edge_thing[current]:
                    if(edges[(current, k)] == str[j]):
                        current = k
                        flag = True
                        break
            if(not flag):
                x = nodes[len(nodes)-1]+1
                nodes.append(x)
                edges[(current, x)] = str[j]
                if(current not in edge_thing):
                    edge_thing[current] = []
                edge_thing[current].append(x)
                current = nodes[len(nodes)-1]
    nodes.append(nodes[len(nodes)-1]+1)
    return nodes, edges, edge_thing

"""
Determine a maximal nonbranching path in a graph
Takes in edges dict and list of vertices
Returns list of all maximal non-branching paths
"""
def maximalNonBranching(edges, vertices):
    paths = []
    indegree = dict()
    outdegree = dict()
    visited = dict()
    for v in vertices:
        indegree[v] = 0
        outdegree[v] = 0
    for key, value in edges.items():
        for j in range(len(value)):
            indegree[value[j]] += 1
            outdegree[key] += 1
    for v in vertices:
        if not (outdegree[v] == 1 and indegree[v] == 1):
            if(outdegree[v] > 0):
                for i in edges[v]:
                    nonbranch = [v]
                    visited[v] = 1
                    while (outdegree[i] == 1 and indegree[i] == 1):
                        visited[i] = 1
                        nonbranch.append(i)
                        i = edges[i][0]
                    nonbranch.append(i)
                    visited[i] = 1
                    paths.append(nonbranch)
    for v in vertices:
        if not v in visited:
            cycle = []
            i = v
            while outdegree[i] == 1 and i not in visited:
                cycle.append(i)
                visited[i] = 1
                i = edges[i][0]
            cycle.append(i)
            paths.append(cycle)
    return paths


"""
Create a suffix tree 
Takes in list of edges and list of paths
Returns modified list of edges
"""
def suffixTree(edges, paths):
    for i in paths:
        start = i[0]
        end = i[len(i)-1]
        s = ""
        for j in range(len(i)-1):
            s+=edges[(i[j], i[j+1])]
            del edges[(i[j], i[j+1])]
        edges[(start, end)] = s
    return edges

def matchOne(i, suffix_list, suffix_array):
    low = 0
    high = len(suffix_list)-1
    first = -1
    last = -1
    while(low <= high):
        mid = math.floor((low+high)/2)
        if (i == suffix_list[mid][0:len(i)]):
            first = mid
            high = mid - 1
        if(i > suffix_list[mid][0:len(i)]):
            low = mid+1
        else:
            high = mid-1
    if(first == -1):
        return -1, -1
    low = 0
    high = len(suffix_list)-1
    while(low <= high):
        mid = math.floor((low+high)/2)
        if(suffix_list[mid][0:len(i)] == i):
            last = mid
            low = mid+1
        elif(suffix_list[mid][0:len(i)] > i):
            high = mid-1
        else:
            low = mid+1
    return first, last


"""
Look for read matches with a set of patterns
Takes in list of patterns and input string
Return all positions where input string matches
"""
def match(patterns, string):
    suffix_array = suffixArray(string)
    suffix_list = []
    positions = {}
    for key, value in suffix_array.items():
        suffix_list.append(key)
    print(len(suffix_list))
    for i in patterns:
        first, last = matchOne(i, suffix_list, suffix_array)
        positions[i] = []
        if(first != -1 and last != -1):
            for j in range(first, last+1):
                positions[i].append(suffix_array[suffix_list[j]])
            positions[i].sort()
    return positions


"""
Determine the amount of different characters in a string
Takes in 2 strings one, two
Returns number of diff characters
"""
def stringDiff(one, two):
    diff = 0
    for i in range(len(one)):
        if(i > len(two)-1 or one[i] != two[i]):
            diff+=1
    return diff

"""
Look for read matches with up to a certain amount of mismatches (d)
Takes in list of patterns, string to match, and d mismatches
Returns list of positions
"""
def matchWithMismatch(patterns, string, d):
    suffix_array = suffixArray(string)
    suffix_list = []
    positions = {}
    for key, value in suffix_array.items():
        suffix_list.append(key)
    for i in patterns:
        k = math.floor(len(i)/(d+1))
        positions[i] = []
        for j in range(0, len(i)-k+1):
            new = i[j:j+k]
            first, last = matchOne(new, suffix_list, suffix_array)
            if(first != -1 and last != -1):
                for l in range(first, last+1):
                    pos = suffix_array[suffix_list[l]]
                    subs_1 = pos-j
                    subs_2 = subs_1+len(i)
                    diff = stringDiff(i, string[subs_1:subs_2])
                    if(diff <= d and subs_1 not in positions[i]):
                        positions[i].append(subs_1)
        positions[i].sort()
    return positions


if __name__ == '__main__':
    f = open("input.txt", "r")
    l = f.readlines()
    string = l[0].strip()
    patterns = l[1].split()
    mismatch = int(l[2])
    x = (matchWithMismatch(patterns, string, mismatch))


    f= open("output.txt", "w")
    for key, value in x.items():
        f.write(key+ ":")
        for i in value:
            f.write(" "+str(i))
        f.write("\n")

    
        