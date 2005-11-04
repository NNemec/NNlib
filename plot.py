############################################################################

import sys,os,colorsys

############################################################################

if os.getenv('NUMERIX') == 'numarray':
    sys.argv += ['--numarray']

import pylab
import matplotlib
#from matplotlib.font_manager import FontProperties

if os.getenv('NUMERIX') == 'numarray':
    sys.argv[-1:] = []

from pylab import (
    plot,subplot,legend,figure,
    gcf,gca,axes,title,
    savefig,
    xlabel,ylabel,xlim,ylim,xticks,yticks,loglog,semilogx,semilogy,
    axhline,axvline,
    matshow,imshow,
    subplots_adjust,
)

############################################################################

def zeroline():
    pylab.axhline(color=(.5,.5,.5))

def savesvg():
    assert sys.argv[0][-3:] == '.py'
    pylab.savefig(sys.argv[0][:-3]+".svg")

def savepng(dpi=50):
    assert sys.argv[0][-3:] == '.py'
    pylab.savefig(sys.argv[0][:-3]+".png",dpi=dpi)

def rainbow(n,N):
    assert 0 <= n
    assert n < N
    return colorsys.hsv_to_rgb(n*.75/N,1,.5)

def complexphasecolor(cpx):
    hs = min(abs(cpx),1.)
    v = atan2(cpx.imag,cpx.real)/(2*pi)+0.5
    return colorsys.hsv_to_rgb(hs,hs,v)

def show():
#    if pylab.rcParams['backend'] == 'FltkAgg':
#        pylab.show(True)
#    else:
	pylab.show()

############################################################################

def figtitle(title):
    pylab.figtext(0.5,0.99,
	title,
	horizontalalignment='center',
	verticalalignment='top',
    )

def set_figsize_inches(w,h):
    pylab.gcf().set_figsize_inches(w,h)


_subplot_grid = [1,1,1]

def subplot_setgrid(rows,cols):
    _subplot_grid[0:3] = [rows,cols,0]

def subplot_start():
    _subplot_grid[2] %= _subplot_grid[0]*_subplot_grid[1]
    _subplot_grid[2] += 1
    pylab.subplot(*_subplot_grid)

############################################################################
