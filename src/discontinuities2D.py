'''
1 Point discontinuities and edges
1.2 Discontinuities in flatland: points and edges
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import pywt
import numpy as np
import utils

import pylab

def matrix():
    S = np.zeros(np.array([100, 100]), dtype=float)
    S[49, 49] = 1
    return S      
    
def matrix2():
    S = np.zeros(np.array([100, 100]), dtype=float)
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            if (i-49)*(i-49) + (j-49)*(j-49) <= 10:
                S[i,j] = 1
    return S

def task_fwt2(S):
   A = pywt.wavedec2(S, "db4", mode=pywt.MODES.ppd, level=4)
   utils.draw_coeffs2(utils.concat_coeffs2(A))
   
def task_fft2(S):
   A = np.fft.fft2(S, s=None, axes=(-2, -1))
   utils.draw_coeffs2(A.real)
   
if __name__ == "__main__":
    S1 = matrix()
    S2 = matrix2()
    
    utils.draw_coeffs2(S1)
    
    #pylab.spy(S1)
    #pylab.spy(S2)
    #utils.draw_coeffs2(S1)
    #utils.draw_coeffs2(S2)
    
    task_fwt2(S1)
    task_fft2(S1)
    task_fwt2(S2)
    task_fft2(S2)