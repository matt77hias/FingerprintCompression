'''
2 Wavelet packets
2.4 Best basis (2D)
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import numpy as np
import pylab

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

###############################################################################
# VISUALIZATIONS
###############################################################################

def visualize_cost_threshold():
    t = [0.1, 0.01, 0.001]
    r = np.arange(-0.2, 0.2, 0.001)
    css = np.zeros((3, r.shape[0]))
    for j in range(len(t)) :
        costf = cost_threshold(t[j])
        for i in range(r.shape[0]):
            css[j,i] = costf(np.array(r[i]))
            
    f, (ax1, ax2, ax3) = pylab.subplots(3, sharex=True, sharey=True)
    ax1.plot(r, css[2], label='Threshold=0.001')
    ax2.plot(r, css[1], label='Threshold=0.01')
    ax3.plot(r, css[0], label='Threshold=0.1')
    f.subplots_adjust(hspace=0)
    pylab.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    pylab.xlabel("Coefficient value")
    pylab.ylabel("Cost")
    pylab.ylim(-1.5,1.5)
    pylab.show()
    
def visualize_cost_shannon():
    r = np.arange(-2, 2, 0.01)
    cs = np.zeros(r.shape)
    costf = cost_shannon
    for i in range(r.shape[0]):
        cs[i] = costf(np.array([r[i]]))
    
    pylab.figure()
    pylab.plot(r, cs)
    pylab.xlabel("Coefficient value")
    pylab.ylabel("Cost")
    pylab.legend(loc=2)
    pylab.show()
    
if __name__ == "__main__": 
    visualize_cost_threshold()
    visualize_cost_shannon()