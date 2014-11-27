'''
1 Point discontinuities and edges
1.1 What's in a discontinuity? Infinitely many frequencies
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import pywt
import numpy as np
import utils

def h(x):
    if x<0:
        return 1+x
    else:
        return -1+x   

def samples(nb_samples= 1000):
    t = (np.linspace(-1,1,nb_samples+1))[:-1]
    return np.vectorize(h)(x=t)
    
def task_fwt(S):
    A = pywt.wavedec(S, "db4", mode=pywt.MODES.ppd, level=4)
    utils.draw_coeffs(utils.concat_coeffs(A))  
    
def task_fft(S):
    A = np.fft.fft(S, n=None, axis=-1)
    utils.draw_coeffs(A)

if __name__ == "__main__":
    S = samples()
    task_fwt(S)
    task_fft(S)