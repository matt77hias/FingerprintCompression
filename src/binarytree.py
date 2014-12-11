'''
2 Wavelet packets
2.4 Best basis (1D)
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import cost
import heapq
import node
import numpy as np
import pywt

###############################################################################
# ANALYSIS ALGORITHM FUNCTIONS
############################################################################### 
       
def wp(S, costf, wavelet="db4", mode=pywt.MODES.ppd, level=None):
    '''
    Returns the 1D discrete wavelet packet transformation, with the best basis according
    to the given cost function, for the given 1D input signal.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param costf:      The (single parameter) cost function that must be used while
                      searching for the best basis.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @param level:     Number of decomposition steps to perform. If the level is None, then the
                      full decomposition up to the level computed with dwt_max_level() function for
                      the given data and wavelet lengths is performed.
    @return:          A list containing the nodes of the 1D discrete wavelet packet transformation,
                      with the best basis according to the given cost function, for the given input signal. 
    '''
    if (level == None):
        level = pywt.dwt_max_level(S.shape[0], pywt.Wavelet(wavelet))
    
    #Data collection step
    Nodes = collect(S, wavelet=wavelet, mode=mode, level=level)
    #Dynamic programming upstream traversal
    mark(Nodes, costf)
    #node.print_nodes(Nodes)
    #Dynamic programming downstream traversal
    Result = []
    traverse(Nodes[0][0], Nodes, Result)
    traverse(Nodes[0][1], Nodes, Result)
    return sorted(Result, cmp=node.compare_low_level_first, reverse=False)
                     
def collect(S, wavelet, mode, level):
    '''
    Returns the full binary tree of wavelet packets.
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
    @return:          The full binary tree of wavelet packets.
    '''
    Nodes = [[] for i in range(level)]
    (Cl, Cr) = pywt.dwt(S, wavelet=wavelet, mode=mode)
    Nodes[0] = [node.Node(Cl, 0, 0), node.Node(Cr, 0, 1)]
    for l in range(0, level-1):
        Parents = Nodes[l]
        Childs = []
        for p in range(len(Parents)):
            (Cl, Cr) = pywt.dwt(Parents[p].C, wavelet=wavelet, mode=mode)
            Childs.append(node.Node(Cl, l+1, 2*p))
            Childs.append(node.Node(Cr, l+1, 2*p+1))
        Nodes[l+1] = Childs 
    return Nodes
    
def mark(Nodes, costf):
    '''
    Marks every node of nodes with the best cost seen so far. 
    @param Nodes:     List containing the nodes of the 1D discrete wavelet packet
                      transformation.
    @param costf:      The (single parameter) cost function that must be used while
                      searching for the best basis.
    '''
    for p in range(len(Nodes[-1])):
        Node = Nodes[-1][p]
        cp = costf(Node.C)
        Node.cost = cp
        Node.best = cp
    for l in range(len(Nodes)-2, -1, -1):
        for p in range(len(Nodes[l])):
            Node = Nodes[l][p]
            cc = Nodes[l+1][2*p].best + Nodes[l+1][2*p+1].best
            cp = costf(Node.C)
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
    @param Nodes:     List containing the nodes of the 1D discrete wavelet packet
                      transformation.
    @param Result:    Buffer containing the nodes traversed so far that belong
                      to the best basis.
    '''
    if (Node.best == Node.cost):
        Result.append(Node)
    else:
        i = Node.level + 1
        j = 2 * Node.index
        traverse(Nodes[i][j], Nodes, Result)
        traverse(Nodes[i][j+1], Nodes, Result)                   

###############################################################################
# SYNTHESIS ALGORITHM FUNCTIONS
############################################################################### 
       
def iwp(Nodes, wavelet="db4", mode=pywt.MODES.ppd):
    '''
    Returns the inverse 1D discrete wavelet packet transformation for the given
    list containing the nodes of the 1D discrete wavelet packet transformation.
    @param Nodes:     List containing the nodes of the 1D discrete wavelet packet
                      transformation.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @return:          The inverse 1D discrete wavelet packet transformation for the given
                      list containing the nodes of the 1D discrete wavelet packet transformation.
    '''
    heapq.heapify(Nodes)
    while len(Nodes) != 1:
        Node1 = heapq.heappop(Nodes)
        Node2 = heapq.heappop(Nodes)
        S = pywt.idwt(Node1.C, Node2.C, wavelet=wavelet, mode=mode, correct_size=True)
        Merged = node.Node(S, (Node1.level-1), (Node1.index / 2))
        heapq.heappush(Nodes, Merged)
    return Nodes[0].C

###############################################################################
# TESTS
############################################################################### 
                                         
def h(x):
    if x<0:
        return 1+x
    else:
        return -1+x   

def samples(nb_samples=1000):
    t = (np.linspace(-1,1,nb_samples+1))[:-1]
    return np.vectorize(h)(x=t)  
      
def get_smooth_func():
    time = np.linspace(1,40,1000)
    freq = 500
    freqs = 8000
    return np.sin(2*np.pi*freq/freqs*time)     
      
if __name__ == "__main__":
    S = get_smooth_func()
    #S = samples()
    Nodes=wp(S, cost.cost_shannon)
    node.print_flattened_nodes(Nodes)
    R = iwp(Nodes)