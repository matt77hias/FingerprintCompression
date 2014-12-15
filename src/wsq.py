import node
import quadtree
import pywt

###############################################################################
# ANALYSIS ALGORITHM FUNCTIONS
###############################################################################  

def subband_decompose(S, wavelet="db4", mode=pywt.MODES.ppd):
    '''
    Returns the 2D subband decomposition for fingerprints for the given 2D input signal.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @return:          A list containing the nodes of the 2D subband decomposition for fingerprints
                      for the given input signal. 
    '''
    #Data collection step
    #Optimization: don't construct the full tree, but lazy initialize
    Nodes = quadtree.collect(S, wavelet=wavelet, mode=mode, level=5)
    #node.print_nodes(Nodes)
    #Downstream traversal
    Result = []
    traverse(Nodes[0][0], Nodes, Result)
    traverse(Nodes[0][1], Nodes, Result)
    traverse(Nodes[0][2], Nodes, Result)
    traverse(Nodes[0][3], Nodes, Result)
    return sorted(Result, cmp=node.compare_low_level_first, reverse=False)
        
def traverse(Node, Nodes, Result):
    '''
    Traverses the given node.
    The node will be added to the result if it belongs to the best basis.
    Otherwise the node childs will be traversed recursively.
    @param Node:      The current node to traverse.
    @param Nodes:     List containing the nodes of the 2D subband decomposition
                      for fingerprints.
    @param Result:    Buffer containing the nodes traversed so far that belong
                      to the best basis.
    '''
    i = Node.level
    j = Node.index
    
    isBottom = False
    if (i==1):
        if(j/4 !=0 or j==3):
            isBottom = True
    elif (i==3 and j!=0):
        isBottom = True
    elif (i == 4):
        isBottom = True
    
    if (isBottom):
        Result.append(Node)
    else: 
        i = i + 1
        j = 4 * j
        traverse(Nodes[i][j], Nodes, Result)
        traverse(Nodes[i][j+1], Nodes, Result)
        traverse(Nodes[i][j+2], Nodes, Result)  
        traverse(Nodes[i][j+3], Nodes, Result)
        
###############################################################################
# SYNTHESIS ALGORITHM FUNCTIONS
###############################################################################
        
def isubband_decompose(Nodes, wavelet="db4", mode=pywt.MODES.ppd):
    '''
    Returns the inverse 2D subband decomposition for fingerprints for the given
    list containing the nodes of the 2D discrete wavelet packet transformation.
    @param Nodes:     List containing the nodes of the 2D subband decomposition
                      for fingerprints
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodic-padding.
    @return:          The inverse 2D subband decomposition for fingerprints for the given
                      list containing the nodes of the 2D discrete wavelet packet transformation.
    '''
    return quadtree.iwp2(Nodes, wavelet=wavelet, mode=mode)

###############################################################################
# TESTS
###############################################################################      

import configuration as c
import cv2
def fingerprint():
    fname = c.get_dir_fingerprints() + "cmp00001.pgm"
    return cv2.imread(fname, 0)   
      
if __name__ == "__main__":
    S = fingerprint()
    Nodes=subband_decompose(S)
    node.print_flattened_nodes(Nodes)