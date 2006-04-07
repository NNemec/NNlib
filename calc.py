import os, sys
from math import *
import numpy
from numpy import *
import scipy.linalg
from scipy.linalg import *
#from numpy.linalg import *

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

def inv(X):
    return X.I # Matrix(numpy.linalg.inv(M))

def eigvals(X):
#    print "eigvals:",
#    X = asarray_chkfinite(X)
#    geev, = scipy.linalg.lapack.get_lapack_funcs(('geev',),(X,))
#    print geev.module_name
    return scipy.linalg.eigvals(X)

def integral(x,y):
    return .5*sum((x[1:]-x[:-1])*(y[1:]+y[:-1]))

#def norm(X):
#    v = asarray(X).flatten()
#    return abs(dot(v,conj(v))**.5)
