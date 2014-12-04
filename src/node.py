'''
2 Wavelet packets
Nodes and node utilities
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

###############################################################################
# PRINT FUNCTIONS
############################################################################### 

def print_nodes(Nodes, fname='nodes.txt'):
    '''
    Prints the given non-flattened nodes to a file with the given file name.
    @param Nodes:     List containing the non-flattened nodes of the (1D or 2D)
                      discrete wavelet packet transformation.
    @param fname:     The filename.
    '''
    f = open(fname,'w')
    for i in range(len(Nodes)):
        for j in range(len(Nodes[i])):
            Node = Nodes[i][j]
            string = "(" + str(Node.level) + "," + str(Node.index) + ")::[" + str(Node.cost) +"|" + str(Node.best) + "]"
            f.write("%s\n" % string) 
    f.close()

def print_flattened_nodes(Nodes, fname='result.txt'):
    '''
    Prints the given flattened nodes to a file with the given file name.
    @param Nodes:     List containing the flattened nodes of the (1D or 2D)
                      discrete wavelet packet transformation.
    @param fname:     The filename.
    '''
    f = open(fname,'w')
    for i in range(len(Nodes)):
        Node = Nodes[i]
        string = "(" + str(Node.level) + "," + str(Node.index) + ")::[" + str(Node.cost) +"|" + str(Node.best) + "]"
        f.write("%s\n" % string) 
    f.close()                                       

###############################################################################
# NODES
###############################################################################  
                                                                                                                                                 
class Node:
    def __init__(self, C, level, index):
        '''
        Creates a new node.
        @param C:         The coefficients for the new node
        @param level:     The level for the new node
        @param index:     The index for the new node
        @note:            level 0 corresponds with the first
                          decomposition in this implementation
                          and not with the original signal as
                          opposed to the assignment.
        '''
        self.C = C
        self.level = level
        self.index = index
        self.cost = -1
        self.best = -1
    
    def __cmp__(self, other):
        '''
        Compares this node against the given node.
        (high levels first, low indices first)
        @param other:     The node to compare this node to.
        '''
        return compare_high_level_first(self, other)

def compare_low_level_first(Node1, Node2) :
    '''
    Compares the first given node against the second
    given node. (low levels first, low indices first)
    @param Node1:     The first node
    @param Node2:     The second node
    '''
    if Node1.level < Node2.level:
        return -1 
    if Node1.level > Node2.level:
        return 1
    if Node1.index < Node2.index:
        return -1
    if Node1.index > Node2.index:
        return 1
    else: 
        return 0
    
def compare_high_level_first(Node1, Node2) :
    '''
    Compares the first given node against the second
    given node. (high levels first, low indices first)
    @param Node1:     The first node
    @param Node2:     The second node
    '''
    if Node1.level > Node2.level:
        return -1 
    if Node1.level < Node2.level:
        return 1
    if Node1.index < Node2.index:
        return -1
    if Node1.index > Node2.index:
        return 1
    else: 
        return 0 