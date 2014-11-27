'''
Utilities
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''

import numpy as np
import pylab

def concat_coeffs(A):
    return reduce(np.append, A[1:], A[0])
    
def concat_coeffs2(A):
    return reduce(combine, A[1:], A[0])
    
def combine(S, (H, V, D)):
    (Sx, Sy) = S.shape
    (Hx, Hy) = H.shape
    (Vx, Vy) = V.shape
    (Dx, Dy) = D.shape
    
    print("----")
    print(S.shape)
    print(H.shape)
    print(V.shape)
    print(D.shape)
    
    NS = np.zeros((Sx+Dx, Sy+Dy))
    print(NS.shape)
    for i in range(Sx):
        for j in range(Sy):
            NS[i,j] = S[i,j]
    for i in range(Hx):
        for j in range(Hy):
            NS[i+Sx,j] = H[i,j]
    for i in range(Vx):
        for j in range(Vy):
            NS[i,j+Sy] = V[i,j]
    for i in range(Dx):
        for j in range(Dy):
            NS[i+Sx,j+Sy] = D[i,j]
    return NS   
    
    
def draw_coeffs(C):
    pylab.figure()
    pylab.semilogy(abs(C), '-')
    pylab.show()
    
def draw_coeffs2(C):
    pylab.figure()
    pylab.imshow(C, interpolation='nearest', cmap=pylab.cm.ocean, extent=(0.5,np.shape(C)[0]+0.5,0.5,np.shape(C)[1]+0.5))
    pylab.colorbar()
    pylab.show()