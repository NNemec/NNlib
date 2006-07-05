#!/usr/bin/env python

from calc import *
from plot import *

def hsv_to_rgb(h, s, v):
#    v=1,s=1,h=0 -> i=0 f=0 p=0 q=1 t=0 -> (1,0,0) red
#	    h=1 -> i=1 (1,1,0)
#    print "h:",h
#    print "s:",s
#    print "v:",v
    f = abs(((h%1.0)*6.0+1.0) % 2.0 - 1.0)
#    print "f:",f

    max_val = v*(1.0)
    mid_val = v*(1.0-s*(1.0-f))
    min_val = v*(1.0-s)

#    print "max_val:",max_val
#    print "mid_val:",mid_val
#    print "min_val:",min_val

    h = asarray(h)
#    print "h:",h
    max_color = (h*3.0 + 0.5).astype(Int) % 3   # red=0, green=1, blue=2
    min_color = (h*3.0 + 2.0).astype(Int) % 3
    mid_color = 3 - (max_color + min_color)
    
#    print "max_color:",max_color
#    print "mid_color:",mid_color
#    print "min_color:",min_color

    rgb = zeros(mid_val.shape+(3,),'d')

    for c in range(3):
	rgb[...,c] = select(
	    [ max_color == c, mid_color == c, min_color == c],
	    [ max_val, mid_val, min_val],
	)

    return rgb[...,0],rgb[...,1],rgb[...,2]

def cpx_to_rgb(cpx):
    Z = abs(cpx)
    re = cpx.real / (Z+1e-5)
    im = cpx.imag / (Z+1e-5)
    
    r =  0.5  * re +      0.0 * im + 0.5
    g = -0.25 * re +  3**.5/4 * im + 0.5
    b = -0.25 * re + -3**.5/4 * im + 0.5
    
    return r,g,b
    

class HSVColormap(matplotlib.colors.Colormap):
    def __init__(self):
	self.N=256
	self.monochrome=False
	self.name='HSV'

    def __call__(self,X,alpha=1.0):
        alpha = min(alpha, 1.0) # alpha must be between 0 and 1
        alpha = max(alpha, 0.0)
        if isinstance(X, (int, float)):
            vtype = 'scalar'
            xa = array([X])
        else:
            vtype = 'array'
            xa = asarray(X)

	r,g,b = hsv_to_rgb(xa,xa,xa)
	
        rgba = zeros(xa.shape+(4,), Float)
        rgba[...,0] = r
        rgba[...,1] = g
        rgba[...,2] = b
        rgba[...,3] = alpha

        if vtype == 'scalar':
            return tuple(rgba[0,:])
        else:
	    return rgba

class SmoothColormap(matplotlib.colors.Colormap):
    def __init__(self,basecolor=(1.0,0.0,0.0)):
	self.N=256
	self.monochrome=False
	self.name='Smooth'
	self.basecolor = array(basecolor)
	assert self.basecolor.shape == (3,)

    def __call__(self,X,alpha=1.0):
        alpha = min(alpha, 1.0) # alpha must be between 0 and 1
        alpha = max(alpha, 0.0)
        if isinstance(X, (int, float)):
            vtype = 'scalar'
            xa = array([X])
        else:
            vtype = 'array'
            xa = asarray(X)

	bc = self.basecolor[:]
	for i in range(len(xa.shape)-1):
	    bc = bc[None,...]

	rgb = (1.0 - (xa[...,None])*(1.0-bc[...,:]))*(1.0 - xa[...,None])

        rgba = zeros(xa.shape+(4,), Float)
        rgba[...,0:3] = rgb
        rgba[...,3] = alpha

        if vtype == 'scalar':
            return tuple(rgba[0,:])
        else:
	    return rgba
	    
    def is_gray(self):
	return False


if __name__ == '__main__':
    x = linspace(0,1,97)
#y = linspace(0,1,2)

    z = x[None,:]# * y[None,:]
#r,g,b = hsv_to_rgb(z,z,z)
#rgb = concatenate((r[:,:,None],g[:,:,None],b[:,:,None],),axis=2)

    imshow(
        z,
    #    extent=(x[0],x[-1],y[0],y[-1]),
        cmap=SmoothColormap(basecolor=(1.0,0.0,1.0)),
        origin='lower',
        interpolation='nearest',
    )

    show()
