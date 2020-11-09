# UPGMA
# Reads a fasta file and prints the sequences in newick format

import math
import numpy as np
import pandas as pd

def getlistofseq(val):
    listofnames = val[0::2]
    listofseq = val[1::2]
    return {'listofnames':listofnames, 'listofseq': listofseq}

def calculatep(seq1, seq2):
    a = b = 0
    for i in range(0, len(seq1)):
        if seq1[i] == ('-') or seq2[i] == ('-'):
            a += 1
        if seq1[i] == seq2[i] == ('-'):
            b += 1
    s= d = 0
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:
            s += 1
        if seq1[i] != seq2[i]:
            d += 1
    S = s - b
    D = d - a
    p = (D / (S + D))
    return p

def calculateKfromp(p):
    return (-3/4)*math.log(1-4*p/3)

def makeemptymatrix(n):
    m = []
    for i in range(n):
        m.append([])
        for j in range(n):
            m[i].append(0)
    assert len(m) == n
    assert len(m[0]) == n
    return m

def pick2(m):
    min_val = float("inf")
    deli, delj = -1, -1

    for i in range(0, len(m)):
        for j in range(0,len(m[i])):
            if (m[i][j] < min_val) and (i != j):
                min_val = m[i][j]
                deli, delj = i, j
    return min_val, deli, delj

def buildnewnode(i, j, dis):
    return "(" + i + "," + j + ':' + str(dis) + ")" 

file = open ("Myotis_aligned.fa",'r')
fasta = file.read()
x = fasta.split()
a = len(x)
l = getlistofseq(x)
n = len(l['listofseq'])
matrix = makeemptymatrix(len(l['listofseq']))

listofnames = []
for i in range(len(l['listofnames'])):
    a = l['listofnames'][i]
    b = a[12:]
    listofnames.append(b)

Klist = []
for i in range(0, n):
    for j in range(0, n):
        seq1 = l['listofseq'][i]
        seq2 = l['listofseq'][j]
        t = [listofnames[i], listofnames[j]]
        if i != j:
            p = calculatep(seq1,seq2)
            K = calculateKfromp(p)
            matrix[i][j] = K
            Klist.append(K)

row_label = listofnames.copy()

while len(row_label) > 1:
    (K,deli,delj) = pick2(matrix)
    
    dis = K/2

    Node1 = row_label[deli]
    Node2 = row_label[delj]
    newnode = buildnewnode(Node1, Node2, dis)
    row_label.append(newnode)


    matrix.append([])
    newrowindex = len(matrix)-1
    for i in range(len(matrix[0])):
        matrix[newrowindex].append(0)
    for i in range(len(matrix)):
        matrix[i].append(0)

    for i in range (0, len(row_label)):
        for j in range(0, len(row_label)):
            if (matrix[i][j] == 0) and (i!= j):
                matrix[i][j] = matrix[j][i] = (matrix[i][deli]+matrix[i][delj])/2

    matrix.remove(matrix[delj])
    matrix.remove(matrix[deli])

    for j in matrix:
        del j[delj]
        del j[deli]
    
    row_label.remove(Node1)
    row_label.remove(Node2)
    
print (row_label[0] + ';')


file = open ("Newickstring_result.txt",'w')
value = (row_label[0] + ';')
file.write(str(value))
file.close()
