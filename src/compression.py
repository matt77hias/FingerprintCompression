'''
3 Fingerprint compression
3.1 Compression
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''
import cost
import numpy as np
import pywt
import quadtree
import utils
import wsq

###############################################################################
# COMPRESSION FUNCTIONS
############################################################################### 

def compress_dwt2(S, fraction, wavelet="db4", mode=pywt.MODES.ppd, level=4, stats=[]):
    '''
    Computes the 2D discrete wavelet transformation for the given 2D input signal.
    Sets all coefficients with an absolute value below the threshold * maximum of the absolute
    values of the coefficients to zero.
    Returns the inverse 2D discrete wavelet transformation for the modified coefficients
    of the 2D discrete wavelet transformation.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param fraction:  The fraction.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodization.
    @param level:     Number of decomposition steps to perform.
    @return:          The inverse 2D discrete wavelet transformation for the modified coefficients
                      of the 2D discrete wavelet transformation.
    '''
    # 2D discrete wavelet transform
    A = pywt.wavedec2(S, wavelet=wavelet, mode=mode, level=level)
    
    # Compression
    maximum = np.amax(abs(A[0]))
    for (CH, CV, CD) in A[1:]:
        maximum = max(maximum, np.amax(abs(CH)), np.amax(abs(CV)), np.amax(abs(CD)))
    threshold = fraction * maximum
    B = [pywt.thresholding.hard(A[0], threshold, 0)]
    for (CH, CV, CD) in A[1:]:
        CCH = pywt.thresholding.hard(CH, threshold, 0)
        CCV = pywt.thresholding.hard(CV, threshold, 0)
        CCD = pywt.thresholding.hard(CD, threshold, 0)
        B.append((CCH, CCV, CCD))
        
    n = utils.number_of_large_coeffs(utils.concat_coeffs2(B), threshold=threshold)
    stats.append(n)
        
    # 2D inverse discrete wavelet transform
    return pywt.waverec2(B, wavelet=wavelet, mode=mode)
    
def compress_wp2(S, fraction, costf=cost.cost_shannon, wavelet="db4", mode=pywt.MODES.ppd, level=4, stats=[]):
    '''
    Computes the 2D discrete wavelet packet transformation, with the best basis according
    to the given cost function, for the given 2D input signal.
    Sets all coefficients with an absolute value below the threshold * maximum of the absolute
    values of the coefficients to zero.
    Returns the inverse 2D discrete wavelet packet transformation for the modified coefficients
    of the 2D discrete wavelet packet transformation.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param fraction:  The fraction.
    @param costf:      The (single parameter) cost function that must be used while
                      searching for the best basis.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodization.
    @param level:     Number of decomposition steps to perform.
    @return:          The inverse 2D discrete wavelet packet transformation for the modified coefficients
                      of the 2D discrete wavelet packet transformation.
    '''
    # 2D discrete wavelet packet transform
    Nodes = quadtree.wp2(S, costf, wavelet=wavelet, mode=mode, level=level)
    
    # Compression
    maximum = -1
    for Node in Nodes:
        maximum = max(maximum, np.amax(abs(Node.C)))
    threshold = fraction * maximum
    for Node in Nodes:
        Node.C = pywt.thresholding.hard(Node.C, threshold, 0)

    n = 0
    for Node in Nodes:   
        n = n + utils.number_of_large_coeffs(Node.C, threshold=threshold)
    stats.append(n)
    
    # 2D inverse discrete wavelet packet transform
    return quadtree.iwp2(Nodes, wavelet=wavelet, mode=mode)
    
def compress_sd(S, fraction, wavelet="db4", mode=pywt.MODES.ppd, stats=[]):
    '''
    Computes the subband decomposition for fingerprints for the given 2D input signal.
    Sets all coefficients with an absolute value below the threshold * maximum of the absolute
    values of the coefficients to zero.
    Returns the inverse subband decomposition for fingerprints for the modified coefficients
    of the subband decomposition for fingerprints.
    @param S:         Input signal.
                      Both single and double precision floating-point data types are supported
                      and the output type depends on the input type. If the input data is not
                      in one of these types it will be converted to the default double precision
                      data format before performing computations.
    @param fraction:  The fraction.
    @param wavelet:   Wavelet to use in the transform. 
                      This must be a name of the wavelet from the wavelist() list.
    @param mode:      Signal extension mode to deal with the border distortion problem.
                      The default mode is periodization.
    @return:          The inverse subband decomposition for fingerprints for the modified coefficients
                      of the subband decomposition for fingerprints.
    '''
    # 2D discrete wavelet packet transform
    Nodes = wsq.sd(S, wavelet=wavelet, mode=mode)
    
    # Compression
    maximum = -1
    for Node in Nodes:
        maximum = max(maximum, np.amax(abs(Node.C)))
    threshold = fraction * maximum
    for Node in Nodes:
        Node.C = pywt.thresholding.hard(Node.C, threshold, 0)
     
    n = 0
    for Node in Nodes:   
        n = n + utils.number_of_large_coeffs(Node.C, threshold=threshold)
    stats.append(n)
    
    # 2D inverse discrete wavelet packet transform
    return wsq.isd(Nodes, wavelet=wavelet, mode=mode)

###############################################################################
# COMPRESSION UTILITIES
############################################################################### 
   
def mse(S1, S2):
    '''
    Returns the mean squared error of the compressed 2D signal S2
    against the original 2D signal S1.
    @param S1:        The original 2D signal
    @param S2:        The compressed 2D signal
    '''
    D = S1-S2
    return (float(np.sum(np.multiply(D, D)))) / (D.shape[0]*D.shape[1])
    
def best_fit(S1, S2):
    (m, n) = S1.shape
    (p, q) = S2.shape
    bi = bj = -1
    best = np.inf
    for i in range(p - m + 1):
        for j in range(q - n + 1):
            error = mse(S1, S2[i:i+m,j:j+n])
            if error < best:
                best = error
                bi = i
                bj = j
    return (S2[bi:bi+m,bj:bj+n], best)

###############################################################################
# TESTS
###############################################################################

import configuration as c
import cv2
import pylab

write_intermediate_results = True

# Note that it would of course be cleaner to change the fraction
# multiple times between the analysis and the synthesis
# but this is just a test method
def compare(fname, fractions, wavelet="db4", mode=pywt.MODES.ppd, level=4):
    stats_dwt2 = []   
    stats_wp2_s = []
    stats_wp2_t = []
    
    S = 255 - cv2.imread(fname, 0)
    E1 = np.zeros(fractions.shape)
    E2 = np.zeros(fractions.shape)
    E3 = np.zeros(fractions.shape)
    i = 0
    for f in fractions:
        R1 = compress_dwt2(S, f, wavelet=wavelet, mode=mode, level=level, stats=stats_dwt2)[level:-level,level:-level]
        R2 = compress_wp2(S, f, costf=cost.cost_shannon, wavelet=wavelet, mode=mode, level=level, stats=stats_wp2_s)[level:-level,level:-level]
        R3 = compress_wp2(S, f, costf=cost.cost_threshold(0.01), wavelet=wavelet, mode=mode, level=level, stats=stats_wp2_t)[level:-level,level:-level]
        R = S[level:-level,level:-level]
        (R1, e1) = best_fit(R, R1)
        (R2, e2) = best_fit(R, R2)
        (R3, e3) = best_fit(R, R3)
        
        if write_intermediate_results:
            S1 = 255 - np.array(R1, dtype=np.uint8)
            S2 = 255 - np.array(R2, dtype=np.uint8)
            S3 = 255 - np.array(R3, dtype=np.uint8)   
            cv2.imwrite(str(i) + "_dwt_" + str(f) + " " + str(e1) + ".png", S1)
            cv2.imwrite(str(i) + "_wp_s_" + str(f) + " " + str(e2) + ".png", S2)
            cv2.imwrite(str(i) + "_wp_t_" + str(f) + " " + str(e3) + ".png", S3)

        E1[i] = e1
        E2[i] = e2
        E3[i] = e3
        i = i + 1
    
    pylab.figure()
    pylab.loglog(fractions, E1, label='DWT')
    pylab.loglog(fractions, E2, label='WP Shannon')
    pylab.loglog(fractions, E3, label='WP Threshold')
    pylab.xlabel("Fraction")
    pylab.ylabel("Mean Square Error")
    pylab.legend(loc=2)
    pylab.show()
    
    pylab.figure()
    pylab.loglog(fractions, stats_dwt2, label='DWT')
    pylab.loglog(fractions, stats_wp2_s, label='WP Shannon')
    pylab.loglog(fractions, stats_wp2_t, label='WP Threshold')
    pylab.xlabel("Fraction")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()
    
def compare2(fname, fractions, costf=cost.cost_shannon, wavelet="db4", mode=pywt.MODES.ppd):  
    stats_sd = []
    stats_wp2_s = []
    stats_wp2_t = []
    
    level = 5
    S = 255 - cv2.imread(fname, 0)
    E1 = np.zeros(fractions.shape)
    E2 = np.zeros(fractions.shape)
    E3 = np.zeros(fractions.shape)
    i = 0
    for f in fractions:
        R1 = compress_sd(S, f, wavelet=wavelet, mode=mode, stats=stats_sd)[level:-level,level:-level]
        R2 = compress_wp2(S, f, costf=cost.cost_shannon, wavelet=wavelet, mode=mode, level=level, stats=stats_wp2_s)[level:-level,level:-level]
        R3 = compress_wp2(S, f, costf=cost.cost_threshold(0.01), wavelet=wavelet, mode=mode, level=level, stats=stats_wp2_t)[level:-level,level:-level]
        R = S[level:-level,level:-level]
        (R1, e1) = best_fit(R, R1)
        (R2, e2) = best_fit(R, R2)
        (R3, e3) = best_fit(R, R3)
        
        if write_intermediate_results:
            S1 = 255 - np.array(R1, dtype=np.uint8)
            S2 = 255 - np.array(R2, dtype=np.uint8)
            S3 = 255 - np.array(R3, dtype=np.uint8)   
            cv2.imwrite(str(i) + "_sd_" + str(f) + " " + str(e1) + ".png", S1)
            cv2.imwrite(str(i) + "_wp_s_" + str(f) + " " + str(e2) + ".png", S2)
            cv2.imwrite(str(i) + "_wp_t_" + str(f) + " " + str(e3) + ".png", S3)

        E1[i] = e1
        E2[i] = e2
        E3[i] = e3
        i = i + 1
    
    pylab.figure()
    pylab.loglog(fractions, E1, label='SD')
    pylab.loglog(fractions, E2, label='WP Shannon')
    pylab.loglog(fractions, E3, label='WP Threshold')
    pylab.xlabel("Fraction")
    pylab.ylabel("Mean Square Error")
    pylab.legend(loc=2)
    pylab.show()
    
    pylab.figure()
    pylab.loglog(fractions, stats_sd, label='SD')
    pylab.loglog(fractions, stats_wp2_s, label='WP Shannon')
    pylab.loglog(fractions, stats_wp2_t, label='WP Threshold')
    pylab.xlabel("Fraction")
    pylab.ylabel("Number of large coefficients")
    pylab.legend(loc=2)
    pylab.show()     
    
if __name__ == "__main__":
    fname = c.get_dir_fingerprints() + "cmp00001.pgm"
    fractions = np.append([0.0], np.power(10, np.arange(-20.0, 0.0, 0.5)))
    #fractions = np.append([0.0], np.power(10, np.arange(-5.0, 0.0, 1.0)))
    compare(fname, fractions)
    
    fname = c.get_dir_fingerprints() + "cmp00002.pgm"
    fractions = np.append([0.0], np.power(10, np.arange(-20.0, 0.0, 0.5)))
    #fractions = np.append([0.0], np.power(10, np.arange(-5.0, 0.0, 1.0)))
    #compare2(fname, fractions)