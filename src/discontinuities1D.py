'''
1 Point discontinuities and edges
1.1 What's in a discontinuity? Infinitely many frequencies
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import pywt
import numpy as np
import pylab
import utils

def h(x):
    if x<0:
        return 1+x
    else:
        return -1+x   

def samples(nb_samples=1000):
    t = (np.linspace(-1,1,nb_samples+1))[:-1]
    return np.vectorize(h)(x=t)
    
def plot_large_coeffs(mode=pywt.MODES.ppd, level=4, threshold=0.1):
    t = np.arange(1, 10000, 100)
    wcounts_haar = np.zeros(t.shape)
    wcounts_db4 =  np.zeros(t.shape)
    wcounts_db8 =  np.zeros(t.shape)
    wcounts_sym4 = np.zeros(t.shape)
    wcounts_sym8 = np.zeros(t.shape)
    wcounts_coif4 =np.zeros(t.shape)
    fcounts =      np.zeros(t.shape)
    
    for i in range(t.shape[0]):
        S = samples(t[i])
        wcounts_haar[i] = utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="haar",  mode=mode, level=level)), threshold)
        wcounts_db4[i] =  utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="db4",   mode=mode, level=level)), threshold)
        wcounts_db8[i] =  utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="db8",   mode=mode, level=level)), threshold)
        wcounts_sym4[i] = utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="sym4",  mode=mode, level=level)), threshold)
        wcounts_sym8[i] = utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="sym8",  mode=mode, level=level)), threshold)
        wcounts_coif4[i] =utils.number_of_large_coeffs(utils.concat_coeffs(pywt.wavedec(S, wavelet="coif4", mode=mode, level=level)), threshold)
        fcounts[i] =      utils.number_of_large_coeffs(np.fft.fft(S, n=None, axis=-1), threshold)
    
    pylab.figure()
    pylab.title("threshold = " + str(threshold))
    pylab.plot(t, wcounts_haar, label='haar')
    pylab.plot(t, wcounts_db4, label='db4')
    pylab.plot(t, wcounts_db8, label='db8')
    pylab.plot(t, wcounts_sym4, label='sym4')
    pylab.plot(t, wcounts_sym8, label='sym8')
    pylab.plot(t, wcounts_coif4, label='coif4')
    pylab.xlabel("Number of samples")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()
    
    pylab.figure()
    pylab.title("threshold = " + str(threshold))
    pylab.plot(t, wcounts_haar, label='haar')
    pylab.plot(t, wcounts_db4, label='db4')
    pylab.plot(t, wcounts_db8, label='db8')
    pylab.plot(t, wcounts_sym4, label='sym4')
    pylab.plot(t, wcounts_sym8, label='sym8')
    pylab.plot(t, wcounts_coif4, label='coif4')
    pylab.plot(t, fcounts, label='DFT')
    pylab.xlabel("Number of samples")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()
    
def task_fwt(S, wavelet="db4", mode=pywt.MODES.ppd, level=4):
    A = pywt.wavedec(S, wavelet=wavelet, mode=mode, level=level)
    utils.draw_coeffs(utils.concat_coeffs(A))  
    
def task_fft(S):
    A = np.fft.fft(S, n=None, axis=-1)
    utils.draw_coeffs(A)

if __name__ == "__main__":
    plot_large_coeffs()
    
    S = samples()
    task_fwt(S)
    task_fft(S)