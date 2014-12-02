'''
3 Fingerprint compression
3.1 Compression
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import configuration as c
import cv2
import numpy as np
import pylab
import pywt
import quadtree

def compress_dwt2(S, fraction, wavelet="db4", mode=pywt.MODES.per, level=4):
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
        
    # 2D inverse discrete wavelet transform
    return pywt.waverec2(B, wavelet=wavelet, mode=mode)
    
def compress_wp2(S, fraction, cost=quadtree.cost_shannon, wavelet="db4", mode=pywt.MODES.per, level=4):
    # 2D discrete wavelet packet transform
    Nodes = quadtree.wp2(S, cost, wavelet=wavelet, mode=mode, level=level)
    
    # Compression
    maximum = -1
    for Node in Nodes:
        maximum = max(maximum, np.amax(abs(Node.C)))
    threshold = fraction * maximum
    for Node in Nodes:
        Node.C = pywt.thresholding.hard(Node.C, threshold, 0)
    
    # 2D inverse discrete wavelet packet transform
    return quadtree.iwp2(Nodes, wavelet=wavelet, mode=mode)
    
def mse(S1, S2):
    count = 0
    n = S1.shape[0]
    m = S1.shape[1]
    for i in range(n):
        for j in range(m):
            count = count + (S1[i,j]-S2[i,j])*(S1[i,j]-S2[i,j])
    return float(count) / (n*m)

if __name__ == "__main__":
    S = cv2.imread(c.get_dir_fingerprints() + "cmp00001.pgm", 0)
    cv2.imshow("Original", S)
    fraction = 0.01
    
    R1 = compress_dwt2(S, fraction)
    R2 = compress_wp2(S, fraction)
    print(mse(S, R1))
    print(mse(S, R2))
    
    S1 = np.array(R1, dtype=np.uint8)
    S2 = np.array(R2, dtype=np.uint8)
    #cv2.imwrite("test.pgm", S2)
    cv2.imshow("Reconstructed DWT", S1) 
    cv2.imshow("Reconstructed WP", S2)     