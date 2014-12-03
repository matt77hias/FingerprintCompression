'''
2 Wavelet packets
2.4 Best basis (2D)
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import numpy as np

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
        for c in C.flatten():
           if (abs(c) > threshold):
                cost = cost + 1
        return cost
    return cost_fixed_threshold
        
def cost_shannon(C):
    '''
    Computes the Shannen entropy of a 2D input signal
    @param C:         Input signal.
    '''
    cost = 0
    for c in C.flatten():
        if (c != 0):
            cost = cost - c*c*np.log2(abs(c))
    return cost