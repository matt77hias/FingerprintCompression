import numpy as np
import quadtree
import pywt

def mrProper(S, costf, wavelet="db4", mode=pywt.MODES.ppd):
    '''
    Returns the 2D discrete mr proper wavelet packet transformation, according to the teachings of mr proper,
    for the given 2D input signal.
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
    @return:          A list containing the nodes of the 2D discrete wavelet packet transformation,
                      with the best basis according to the given cost function, for the given input signal. 
    '''
    level=5
    #Data collection step
    Nodes = quadtree.collect(S, wavelet=wavelet, mode=mode, level=level)
    
    #Its going down, im yelling timbeeeeeer
     
def letsAllHaveANiceTime(Node, Nodes, Result):
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
    i = Node.level + 1
    j = 4 * Node.index
    
    isBottom = False
    
    if (i==2):
        if(np.mod(j,4)==3):
            isBottom = True
    elif (i==4 & j!=0):
        isBottom = True
    elif (i==5):
        isBottom = True
 
    if (isBottom):
        Result.append(Node)
    else: 
        letsAllHaveANiceTime(Nodes[i][j], Nodes, Result)
        letsAllHaveANiceTime(Nodes[i][j+1], Nodes, Result)
        letsAllHaveANiceTime(Nodes[i][j+2], Nodes, Result)  
        letsAllHaveANiceTime(Nodes[i][j+3], Nodes, Result)