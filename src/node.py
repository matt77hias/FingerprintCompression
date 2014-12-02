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
    f = open(fname,'w')
    for i in range(len(Nodes)):
        for j in range(len(Nodes[i])):
            Node = Nodes[i][j]
            string = "(" + str(Node.level) + "," + str(Node.index) + ")::[" + str(Node.cost) +"|" + str(Node.best) + "]"
            f.write("%s\n" % string) 
    f.close()

def print_flattened_nodes(Nodes, fname='result.txt'):
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
        self.C = C
        self.level = level
        self.index = index
        self.cost = -1
        self.best = -1

def compare_low_level_first(Node1, Node2) :
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