import os, sys
import numpy
from numpy import *
from numpy.linalg import *

Matrix = numpy.matrix

#import scipy
#from scipy import *
#from scipy.linalg import *
from time import time

#import numarray
#import Numeric
#if os.getenv('NUMERIX') == 'numarray':
#    from numarray.matrix import Matrix
#else:
#    from scipy.linalg import *
#    from Matrix import Matrix

def adj(X):
    return X.H # conj(transpose(X))

def inv(M):
    return X.I # Matrix(numpy.linalg.inv(M))

def eigvals(M):
    return numpy.linalg.eigvals(M)

def integral(x,y):
    return .5*sum((x[1:]-x[:-1])*(y[1:]+y[:-1]))
