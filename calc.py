import os, sys
import scipy
from scipy import *
from scipy.linalg import *

#import numarray
#import Numeric
#if os.getenv('NUMERIX') == 'numarray':
#    from numarray.matrix import Matrix
#else:
#    from scipy.linalg import *
#    from Matrix import Matrix

def adj(X):
    return conj(transpose(X))

def inv(M):
    return Matrix(scipy.linalg.inv(M))

def eigvals(M):
    return scipy.linalg.eigvals(M)

def integral(x,y):
    return .5*sum((x[1:]-x[:-1])*(y[1:]+y[:-1]))
