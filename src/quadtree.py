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

def cost_threshold(threshold):
    '''
    Returns a cost function for computing the number of entries
    of a 1D input signal higher (in absolute value) than the given threshold.
    @param threshold:     The threshold value.
    '''
    def cost_fixed_threshold(C):
        '''
        Computes the number of entries of a 2D input signal
        higher (in absolute value) than the threshold.
        @param C:         Input signal.
        '''
        cost = 0
        for i in range(C.shape[0]):
            for j in range(C.shape[1]):
                if (abs(C[i,j]) > threshold):
                    cost = cost + 1
        return cost
    return cost_fixed_threshold
        
def cost_shannon(C):
    '''
    Computes the Shannen entropy of a 2D input signal
    @param C:         Input signal.
    '''
    cost = 0
    for i in range(C.shape[0]):
         for j in range(C.shape[1]):
            if (C[i,j] != 0):
                cost = cost - C[i,j]*C[i,j]*np.log2(abs(C[i,j]))
    return cost

###############################################################################
# ANALYSIS ALGORITHM FUNCTIONS
###############################################################################        

def wp2(S, cost, wavelet="db4", mode=pywt.MODES.ppd, level=2):
    '''
    Returns the 2D discrete wavelet packet transformation, with the best basis according
    to the given cost function, for the given 2D input signal.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param cost:      The (single parameter) cost function that must be used while
                      searching for the best basis.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @param level:     Number of decomposition steps to perform.
    @return:          A list containing the nodes of the 2D discrete wavelet packet transformation,
                      with the best basis according to the given cost function, for the given input signal. 
    '''
    #Data collection step
    Nodes = collect2(S, wavelet=wavelet, mode=mode, level=level)
    #Dynamic programming upstream traversal
    mark(Nodes, cost)
    node.print_nodes(Nodes)
    #Dynamic programming downstream traversal
    Result = []
    traverse(Nodes[0][0], Nodes, Result)
    traverse(Nodes[0][1], Nodes, Result)
    traverse(Nodes[0][2], Nodes, Result)
    traverse(Nodes[0][3], Nodes, Result)
    return sorted(Result, cmp=node.compare_low_level_first, reverse=False)
                     
def collect2(S, wavelet, mode, level):
    '''
    Returns the full quad tree of wavelet packets.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
    @param level:     Number of decomposition steps to perform. If the level is None, then the
                      full decomposition up to the level computed with dwt_max_level() function for
                      the given data and wavelet lengths is performed.
    @return:          The full quad tree of wavelet packets.
    '''
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
    '''
    Marks every node of nodes with the best cost seen so far. 
    @param Nodes:     List containing the nodes of the 2D discrete wavelet packet
                      transformation.
    @param cost:      The (single parameter) cost function that must be used while
                      searching for the best basis.
    '''
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
    '''
    Traverses the given node.
    The node will be aadded to the result if it belongs to the best basis.
    Otherwise the node childs will be traversed recursively.
    @param Node:      The current node to traverse.
    @param Nodes:     List containing the nodes of the 2D discrete wavelet packet
                      transformation.
    @param Result:    Buffer containing the nodes traversed so far that belong
                      to the best basis.
    '''
    if (Node.best == Node.cost):
        Result.append(Node)
    else:
        i = Node.level + 1
        j = 4 * Node.index
        traverse(Nodes[i][j], Nodes, Result)
        traverse(Nodes[i][j+1], Nodes, Result)
        traverse(Nodes[i][j+2], Nodes, Result)  
        traverse(Nodes[i][j+3], Nodes, Result) 
        
###############################################################################
# SYNTHESIS ALGORITHM FUNCTIONS
###############################################################################        
def iwp2(Nodes, wavelet="db4", mode=pywt.MODES.ppd):
    '''
    Returns the inverse 2D discrete wavelet packet transformation for the given
    list containing the nodes of the 2D discrete wavelet packet transformation.
    @param Nodes:     List containing the nodes of the 2D discrete wavelet packet
                      transformation.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @return:          The inverse 2D discrete wavelet packet transformation for the given
                      list containing the nodes of the 2D discrete wavelet packet transformation.
    '''
    while len(Nodes) != 1:
        Nodes = sorted(Nodes, cmp=node.compare_high_level_first, reverse=False)
        Node1 = Nodes[0]
        Node2 = Nodes[1]
        Node3 = Nodes[2] 
        Node4 = Nodes[3]  
        S = pywt.idwt2((Node1.C, (Node2.C, Node3.C, Node4.C)), wavelet=wavelet, mode=mode)
        Merged = node.Node(S, (Node1.level-1), (Node1.index / 4))
        Nodes = Nodes[4:]
        Nodes.append(Merged)
    return Nodes[0].C
        
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
    Nodes=wp2(S, cost_shannon)
    node.print_flattened_nodes(Nodes)
    R = iwp2(Nodes)