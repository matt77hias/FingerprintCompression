'''
2 Wavelet packets
2.4 Best basis (2D)
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import node
import numpy as np
import pywt

###############################################################################
# COST FUNCTIONS
############################################################################### 

def cost_threshold(C, threshold):
    def cost_fixed_threshold(C):
        cost = 0
        for i in range(C.shape[0]):
            for j in range(C.shape[1]):
                if (abs(C[i,j]) > threshold):
                    cost = cost + 1
        return cost
    return cost_fixed_threshold
        
def cost_shannon(C):
    cost = 0
    for i in range(C.shape[0]):
         for j in range(C.shape[1]):
            if (C[i,j] != 0):
                cost = cost - C[i,j]*C[i,j]*np.log2(abs(C[i,j]))
    return cost

###############################################################################
# ALGORITHM FUNCTIONS
###############################################################################        

def wptree(S, cost, wavelet="db4", mode=pywt.MODES.ppd, level=2):
    #Data collection step
    Nodes = collect(S, wavelet=wavelet, mode=mode, level=level)
    #Dynamic programming upstream traversal
    mark(Nodes, cost)
    node.print_nodes(Nodes)
    #Dynamic programming downstream traversal
    Result = []
    traverse(Nodes[0][0], Nodes, Result)
    traverse(Nodes[0][1], Nodes, Result)
    return sorted(Result, cmp=node.compare, reverse=False)
                     
def collect(S, wavelet, mode, level):
    Nodes = [[] for i in range(level)]
    (CA, (CH, CV, CD)) = pywt.dwt2(S, wavelet=wavelet, mode=mode)
    Nodes[0] = [node.Node(CA, 0, 0), node.Node(CH, 0, 1), node.Node(CV, 0, 2), node.Node(CD, 0, 3)]
    for l in range(0, level-1):
        Parents = Nodes[l]
        Childs = []
        for p in range(len(Parents)):
            (CA, (CH, CV, CD)) = pywt.dwt2(Parents[p].C, wavelet=wavelet, mode=mode)
            Childs.append(node.Node(CA, l+1, 4*p))
            Childs.append(node.Node(CH, l+1, 4*p+1))
            Childs.append(node.Node(CV, l+1, 4*p+2))
            Childs.append(node.Node(CD, l+1, 4*p+3))
        Nodes[l+1] = Childs 
    return Nodes
    
def mark(Nodes, cost):
    for p in range(len(Nodes[-1])):
        Node = Nodes[-1][p]
        cp = cost(Node.C)
        Node.cost = cp
        Node.best = cp
    for l in range(len(Nodes)-2, -1, -1):
        for p in range(len(Nodes[l])):
            Node = Nodes[l][p]
            cc = Nodes[l+1][4*p].best + Nodes[l+1][4*p+1].best + Nodes[l+1][4*p+2].best + Nodes[l+1][4*p+3].best
            cp = cost(Node.C)
            Node.cost = cp
            if cp <= cc:
                Node.best = cp
            else:
                Node.best = cc 
          
def traverse(Node, Nodes, Result):
    if (Node.best == Node.cost):
        Result.append(Node)
        return
    else:
        i = Node.level + 1
        j = 4 * Node.index
        traverse(Nodes[i][j], Nodes, Result)
        traverse(Nodes[i][j+1], Nodes, Result)
        traverse(Nodes[i][j+2], Nodes, Result)  
        traverse(Nodes[i][j+3], Nodes, Result) 
        
###############################################################################
# TESTS
###############################################################################                                          
def matrix(size = 100):
    S = np.zeros(np.array([size, size]), dtype=float)
    half = int(size / 2) - 1
    S[half, half] = 1
    return S      
    
def matrix2(size = 100):
    S = np.zeros(np.array([size, size]), dtype=float)
    half = int(size / 2) - 1
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            if (i-half)*(i-half) + (j-half)*(j-half) <= 10:
                S[i,j] = 1
    return S
      
if __name__ == "__main__":
    S = matrix(64)
    Nodes=wptree(S, cost_shannon)
    node.print_flattened_nodes(Nodes)