import random

"""
Functions to implement genome assembly from many reads using deBruijn graphs
"""

"""
Create de Bruijn graph of adjacent k-mers
Takes in list of kmers, returns dict of edges and vertex list
"""
def deBruijn(kmers):
    vertices = []
    k = len(kmers[0])-1
    d = dict()
    for i in kmers:
        d[i[0:k]] = []
    for i in kmers:
        d[i[0:k]].append(i[1:k+1])
        if(i[0:k] not in vertices):
            vertices.append(i[0:k])
        if(i[1:k+1] not in vertices):
            vertices.append(i[1:k+1])
    d = dict(sorted(d.items()))
    return d, vertices

"""
Create a random walk until a cycle is found through a graph
Takes in edges matrix, a visited list, a start vertex, 
Returns list with path
"""
def randomWalk(edges, start, visited, edgecount):
    path = []
    path.append(start)
    next = random.choice(edges[start])
    while(visited[start][next] == edgecount[start][next]):
        next = random.choice(edges[start])
    path.append(next)
    visited[start][next] += 1
    print("START: ", start, "NEXT: ", next)
    while next != start:
        pot = random.choice(edges[next])
        print("START: ", next, "NEXT: ", pot)
        try:
            while(visited[next][pot] == edgecount[next][pot]):
                pot = random.choice(edges[next])
            visited[next][pot] += 1
            next = pot
            path.append(next)
        except Exception:
            print("FAILED: ", next, pot)
            return
    print(path)
    return path

def checkFilled(edges, visited, num, edgecount):
    for i in range(num):
        if i in edges:
            for j in range(len(edges[i])):
                if visited[i][edges[i][j]] != edgecount[i][edges[i][j]]:
                    return False
    return True


def checkOneVertex(edges, visited, vertex, edgecount):
    print("CHECKING: ", vertex)
    for j in range(len(edges[vertex])):
        if visited[vertex][edges[vertex][j]] != edgecount[vertex][edges[vertex][j]]:
            print("CAN START ON: ", vertex)
            return True
    return False


"""
Look for an Eulerian cycle that visits every edge exactly once
Takes in dict of edges, number of vertices, start node, visited list, and 2D list of edge counts
Returns Eulerian cycle as list
"""
def eulerianCycle(edges, num, start, edgecount):
    visited = [[0 for i in range(num)] for j in range(num)]
    path = randomWalk(edges, start, visited, edgecount)
    print("RETURNED")
    while(not checkFilled(edges, visited, num, edgecount)):
        i = 1
        while(not checkOneVertex(edges, visited, path[i], edgecount)):
            i+=1
        newPath = randomWalk(edges, path[i], visited, edgecount)
        path = path[0:i] + newPath + path[i+1:len(path)]
        print("PATH", path)
    print("RETURNING")
    return path
        

"""
Look for an Eulerian path that visits every edge exactly once
Takes in edges as dict, number of vertices, and 2D list of edge counts
Returns Eulerian path as list
"""
def eulerianPath(edges, num, edgecount):
    indegree = [0 for i in range(num)]
    outdegree = [0 for i in range(num)]
    for i in range(num):
        if(i in edges):
            for j in range(len(edges[i])):
                print(edges[i][j])
                indegree[edges[i][j]]+=1
                outdegree[i]+=1
    start = 0
    end = 0
    for i in range(num):
        if(indegree[i] == outdegree[i]-1):
            start = i
        if(outdegree[i] == indegree[i]-1):
            end = i
    print(start, end)
    
    print("ADDING EDGE: ", end, start)
    if end in edges:
        edges[end].append(start)
    else:
        edges[end] = [start]
    edgecount[end][start]+=1
    print(edges)
    return eulerianCycle(edges, num, start, edgecount)


def reassemble(kmers):
    k = len(kmers[0])-1
    d = dict()
    kmer_to_number = dict()
    number_to_kmer = dict()
    j = -1
    for i in kmers:
        if(not i[0:k] in kmer_to_number):
            print(i[0:k])
            j+=1
            kmer_to_number[i[0:k]] = j
            number_to_kmer[j] = i[0:k]
        if(not i[1:k+1] in kmer_to_number):
            print(i[1:k+1])
            j+=1
            kmer_to_number[i[1:k+1]] = j
            number_to_kmer[j] = i[1:k+1]
        d[kmer_to_number[i[0:k]]] = []

    visited = [[0 for i in range(j+1)] for k in range(j+1)]
    for i in kmers:
        print(kmer_to_number[i[0:k]], kmer_to_number[i[1:k+1]])
        d[kmer_to_number[i[0:k]]].append(kmer_to_number[i[1:k+1]])
        visited[kmer_to_number[i[0:k]]][kmer_to_number[i[1:k+1]]]+=1
    path = eulerianPath(d, j+1, visited)
    i = 2
    string = ""
    string+=(number_to_kmer[path[0]]+number_to_kmer[path[1]][k-1:k])
    print("STRING", string)
    while i < len(path)-1:
        string+=number_to_kmer[path[i]][k-1:k]
        i+=1
    return string

"""
All maximal non branching path in a graph (no intermediate nodes with non-1 indegree and outdegree)
Takes in dict of edges and list of vertices
Returns paths as list of lists
"""  
def maximalNonBranching(edges, vertices):
    paths = []
    indegree = dict()
    outdegree = dict()
    visited = dict()
    print(vertices)
    print(edges)
    for v in vertices:
        indegree[v] = 0
        outdegree[v] = 0
    for key, value in edges.items():
        for j in range(len(value)):
            indegree[value[j]] += 1
            outdegree[key] += 1
    print(outdegree)
    print(indegree)
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
    print(paths)
    for v in vertices:
        if not v in visited:
            cycle = []
            i = v
            while outdegree[i] == 1 and i not in visited:
                cycle.append(i)
                visited[i] = 1
                i = edges[i][0]
            cycle.append(i)
            print(cycle)
            paths.append(cycle)
    return paths

"""
Put it all together
"""
def maximalNonBranchingGenome(lines):
    edges, vertices = deBruijn(lines)
    maximalNon = maximalNonBranching(edges, vertices)
    paths = []
    for i in maximalNon:
        path = ""
        path+=(i[0])
        for j in range(1, len(i)):
            path+=i[j][len(i[j])-1]
        paths.append(path)
    return paths

