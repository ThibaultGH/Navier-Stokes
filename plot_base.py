import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


def read(filename):

    
    source = open(filename, 'r')

    line = source.readlines()

    for i,b in enumerate(line):
        line[i] = b.replace("\n","")

    dim = line[0].split(" ")

    NbNodes = int(dim[0])
    #print(NbNodes)
    NbEle = int(dim[1])
    #print(NbEle)
    
    tab_nodes = np.zeros((NbNodes,2))
    tab_ele = np.zeros((NbEle,3))

    for i in range(0,NbNodes):
        caractere = line[i+1].split(" ")
        #print(caractere)
        for k in range(2):
            tab_nodes[i-1,k] = float(caractere[k])

    for j in range(NbNodes+1,NbEle+NbNodes+1):
        caractere = line[j].split(" ")
        for k in range(3):
            tab_ele[j-NbNodes-1,k] = int(caractere[k])

    return [tab_nodes,tab_ele]
