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

def matrix(size = 100):
    S = np.zeros(np.array([size, size]), dtype=float)
    half = int(size / 2) - 1
    S[half, half] = 1
    return S      
    
def matrix2(size = 100):
    S = np.zeros(np.array([size, size]), dtype=float)
    half = int(size / 2) - 1
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            if (i-half)*(i-half) + (j-half)*(j-half) <= 10:
                S[i,j] = 1
    return S

def plot_large_coeffs2_matrix(mode=pywt.MODES.ppd, level=4, threshold=0.1):
    t = np.arange(1, 200, 10)
    wcounts_haar = np.zeros(t.shape)
    wcounts_db4 =  np.zeros(t.shape)
    wcounts_db8 =  np.zeros(t.shape)
    wcounts_sym4 = np.zeros(t.shape)
    wcounts_sym8 = np.zeros(t.shape)
    wcounts_coif4 =np.zeros(t.shape)
    fcounts =      np.zeros(t.shape)
    
    for i in range(t.shape[0]):
        S = matrix(t[i])
        wcounts_haar[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="haar",  mode=mode, level=level)), threshold)
        wcounts_db4[i] =  utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="db4",   mode=mode, level=level)), threshold)
        wcounts_db8[i] =  utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="db8",   mode=mode, level=level)), threshold)
        wcounts_sym4[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="sym4",  mode=mode, level=level)), threshold)
        wcounts_sym8[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="sym8",  mode=mode, level=level)), threshold)
        wcounts_coif4[i] =utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="coif4", mode=mode, level=level)), threshold)
        fcounts[i] =      utils.number_of_large_coeffs(np.fft.fft2(S, s=None, axes=(-2, -1)).real, threshold)
    
    pylab.figure()
    pylab.title("threshold = " + str(threshold))
    pylab.plot(t, wcounts_haar, label='haar')
    pylab.plot(t, wcounts_db4, label='db4')
    pylab.plot(t, wcounts_db8, label='db8')
    pylab.plot(t, wcounts_sym4, label='sym4')
    pylab.plot(t, wcounts_sym8, label='sym8')
    pylab.plot(t, wcounts_coif4, label='coif4')
    pylab.xlabel("Matrix size")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=1)
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
    pylab.xlabel("Matrix size")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()

def plot_large_coeffs2_matrix2(mode=pywt.MODES.ppd, level=4, threshold=0.1):
    t = np.arange(1, 200, 10)
    wcounts_haar = np.zeros(t.shape)
    wcounts_db4 =  np.zeros(t.shape)
    wcounts_db8 =  np.zeros(t.shape)
    wcounts_sym4 = np.zeros(t.shape)
    wcounts_sym8 = np.zeros(t.shape)
    wcounts_coif4 =np.zeros(t.shape)
    fcounts =      np.zeros(t.shape)
    
    for i in range(t.shape[0]):
        S = matrix2(t[i])
        wcounts_haar[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="haar",  mode=mode, level=level)), threshold)
        wcounts_db4[i] =  utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="db4",   mode=mode, level=level)), threshold)
        wcounts_db8[i] =  utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="db8",   mode=mode, level=level)), threshold)
        wcounts_sym4[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="sym4",  mode=mode, level=level)), threshold)
        wcounts_sym8[i] = utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="sym8",  mode=mode, level=level)), threshold)
        wcounts_coif4[i] =utils.number_of_large_coeffs(utils.concat_coeffs2(pywt.wavedec2(S, wavelet="coif4", mode=mode, level=level)), threshold)
        fcounts[i] =      utils.number_of_large_coeffs(np.fft.fft2(S, s=None, axes=(-2, -1)).real, threshold)
    
    pylab.figure()
    pylab.title("threshold = " + str(threshold))
    pylab.plot(t, wcounts_haar, label='haar')
    pylab.plot(t, wcounts_db4, label='db4')
    pylab.plot(t, wcounts_db8, label='db8')
    pylab.plot(t, wcounts_sym4, label='sym4')
    pylab.plot(t, wcounts_sym8, label='sym8')
    pylab.plot(t, wcounts_coif4, label='coif4')
    pylab.xlabel("Matrix size")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=1)
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
    pylab.xlabel("Matrix size")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()


def task_fwt2(S):
   A = pywt.wavedec2(S, "db4", mode=pywt.MODES.ppd, level=4)
   utils.draw_coeffs2(abs(utils.concat_coeffs2(A)))
   
def task_fft2(S):
   A = np.fft.fft2(S, s=None, axes=(-2, -1))
   utils.draw_coeffs2(abs(A.real))
   
if __name__ == "__main__":
    plot_large_coeffs2_matrix()
    plot_large_coeffs2_matrix2()
    
    S1 = matrix()
    S2 = matrix2()
    utils.draw_coeffs2(S1)
    utils.draw_coeffs2(S2)
    #pylab.spy(S1)
    #pylab.spy(S2)
    task_fwt2(S1)
    task_fft2(S1)
    task_fwt2(S2)
    task_fft2(S2)