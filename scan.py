from numpy import *
import sys
from time import time

class scan_adaptive:
    def __init__(self,
        f,
        initial_xgrid=[],
        period=None,
        xlims=None,
        ylims=None,
        xminstep=None,
        yminstep=None,
        verbose=False,
        sameyscale=False, # use the same y-scale for all values
        masksensitive=None,
        precision=1e-3,
        randomize=0.0,
    ):
        self.verbose = verbose
        self.f = f
        self.period = period
        self.xlims = xlims
        self.ylims = ylims
        self.xminstep = xminstep
        self.yminstep = yminstep
        self.sameyscale = sameyscale
        self.masksensitive = masksensitive
        self.precision = precision
        self.randomize = randomize
        self.x = array([],'d')
        if len(initial_xgrid):
            self.addpoints(initial_xgrid)

    def _initshape(self): # initialization after first value has been red
        yshape = self.y.shape[1:]
        if self.masksensitive is None:
            self.masksensitive = ones(yshape,bool)
        else:
            self.masksensitive = asarray(self.masksensitive).ravel()
            assert masksensitive.shape == yshape
        if self.ylims is not None:
            self.ylims = (asarray(self.ylims[0]).ravel(),asarray(self.ylims[1]).ravel())
            assert self.ylims[0].shape in [(1,),yshape]
            assert self.ylims[1].shape in [(1,),yshape]
        if self.yminstep is not None:
            self.yminstep = asarray(self.yminstep).ravel()
            assert self.yminstep.shape in [(1,),yshape]

    def set_xlims(self,xlims):
        self.xlims = (xlims[0], xlims[1])

    def set_ylims(self,ylims):
        self.ylims = (asarray(ylims[0]).ravel(), asarray(ylims[1]).ravel())
        assert self.ylims[0].shape in [(1,),self.y.shape[1:]]
        assert self.ylims[1].shape in [(1,),self.y.shape[1:]]

    def set_precision(self,precision):
        self.precision = precision

    def debugout(self,str):
        if self.verbose:
            print str,
            sys.stdout.flush()

    def store(self,h5group,gname):
	import tables as pytables
	
        if hasattr(h5group,gname):
            getattr(h5group,gname)._f_remove(recursive=True)
        g = pytables.Group(h5group,gname,new=True)
        filters = pytables.Filters(complevel=9, complib='zlib')
        xatom = pytables.Atom.from_dtype(self.x.dtype)
        yatom = pytables.Atom.from_dtype(self.y.dtype)
        xdata = pytables.CArray(g,'x',shape=self.x.shape,atom=xatom,filters=filters)
        xdata[:] = self.x[:]
        if self.y.shape[1] == 1:
            ydata = pytables.CArray(g,'y',shape=self.y.shape[:1],atom=yatom,filters=filters)
            ydata[:] = self.y[:,0]
        else:
            ydata = pytables.CArray(g,'y',shape=self.y.shape,atom=yatom,filters=filters)
            ydata[:,:] = self.y[:,:]
        h5group._v_file.flush()


    def retrieve(self,h5group,gname,optional=False):
	import tables as pytables
	
        if not hasattr(h5group,gname):
            assert optional
            return
        self.x = getattr(h5group,gname).x[:]
        self.y = getattr(h5group,gname).y[:,...]
        if len(self.y.shape) == 1:
            self.y = self.y[:,None]
        self._initshape()


    def do_scan(self):
        oldN = 0
        while len(self.x) > oldN:
            oldN = len(self.x)
            self.refine_visible()
	nans = isnan(self.y).any(axis=1)
	if any(nans):
	    self.x = self.x[~nans]
	    self.y = self.y[~nans,:]


    def addpoints(self,newxpoints,xlims=None,xminstep=None,internalcall=False):
        if xlims is None:
            xlims = self.xlims
        if xminstep is None:
            xminstep = self.xminstep

        newxpoints = asarray(newxpoints,'d')
        assert len(newxpoints.shape) == 1
        if len(newxpoints) == 0:
            self.debugout("add no points.\n")
            return
        if not internalcall:
            self.debugout("addpoints: old=%i, new=%i, "%(len(self.x),len(newxpoints)))

        if self.period is not None:
            newxpoints = newxpoints % self.period
            if xlims is None:
                xlims = (0.,self.period)
            else:
                xlims = asarray(xlims) % self.period
        if xlims is not None:
            oldcount = len(newxpoints)
            if xlims[0] >= xlims[1] and self.period is not None:
                newxpoints = newxpoints[(newxpoints >= xlims[0]) | (newxpoints <= xlims[1])]
            else:
                newxpoints = newxpoints[(newxpoints >= xlims[0]) & (newxpoints <= xlims[1])]
            if len(newxpoints) < oldcount:
                self.debugout("xlims-cutoff: %i, "%(oldcount - len(newxpoints)))

        if len(self.x) == 0:
            self.x = newxpoints[:1]
            newxpoints = newxpoints[1:]
            y0 = asarray(self.f(self.x[0])).ravel()
            self.y = resize(y0,self.x.shape + y0.shape)
            self.y[0,:] = y0
            self._initshape()

        newxpoints.sort()

        if xminstep is None:
            if xlims is not None:
                xminstep = (xlims[1] - xlims[0]) * self.precision
            else:
                xminstep = (max(self.x[-1],newxpoints[-1]) - min(self.x[0],newxpoints[0])) * self.precision

        slot = self.x.searchsorted(newxpoints)
        oldcount = len(newxpoints)
        newxpoints = newxpoints[
            ((slot == len(self.x)) | (self.x[slot%len(self.x)] - newxpoints >= xminstep))
            &
            ((slot == 0)           | (newxpoints - self.x[slot-1] >= xminstep))
        ]
        if len(newxpoints) < oldcount:
            self.debugout("xminstep-cutoff1: %i, "%(oldcount - len(newxpoints)))

        xdiff = newxpoints[1:]-newxpoints[:-1]
        oldcount = len(newxpoints)
        newxpoints = concatenate((newxpoints[:1],newxpoints[1:][xdiff >= xminstep]))
        if len(newxpoints) < oldcount:
            self.debugout("xminstep-cutoff2: %i, "%(oldcount - len(newxpoints)))

        oldcount = len(newxpoints)

        allx = concatenate((self.x,newxpoints))

        assert allx.dtype.char == 'd'
        allxrounded = (allx/xminstep).round()
        Nduplicate = len(newxpoints) - len(unique(allxrounded[len(self.x):]))
        if Nduplicate > 0:
            self.debugout("duplicate: %i, "%(Nduplicate))
        allxtosort = allxrounded + linspace(0,0.3,len(allxrounded))
        idxsorted = allxtosort.argsort()
        allxsorted = allxrounded[idxsorted]
        idxisunique = concatenate(([True],(allxsorted[1:] != allxsorted[:-1])))
        idxunique = idxsorted[idxisunique | (idxsorted<len(self.x))]
        idxisnew = idxunique >= len(self.x)

        self.x = allx[idxunique]
        self.y = resize(self.y,(len(allx),)+self.y.shape[1:])[idxunique]
        doidx, = nonzero(idxisnew)
        Nrounding = oldcount - len(doidx) - Nduplicate
        if Nrounding > 0:
            self.debugout("rounding: %i, "%(Nrounding))
        if len(doidx) == 0:
            self.debugout("do=0, nothing to do.\n")
        else:
            self.debugout("do=%i ... "%len(doidx))
            starttime = time()
            for i in doidx:
                self.y[i,...] = asarray(self.f(self.x[i])).ravel()
            self.debugout("done. (%g s/point)\n"%((time()-starttime)/len(doidx)))


    def divide(self,
        N=2,
        randomize=False, # ToDo
    ):
        newxpoint = []
        for n in range(1,N):
            newxpoints += list(n*self.x[:-1]/N + (N-n)*self.x[:-1]/N)
        if self.period is not None:
            newxpoints += list(linspace(self.x[-1],self.x[0]+self.period,N+1)[1:-1])
        self.addpoints(newxpoints)


    def find_extrema(self):
        y = self.y[:,self.masksensitive]
        y = concatenate((y[-1:,:],y,y[:1,:]),axis=0)
        diff = y[1:,:] - y[:-1,:]
        extrema = (diff[1:,:] * diff[:-1,:]) <= 0
        if self.period is None:
            extrema[0,:] = False
            extrema[-1,:] = False
        return extrema


    def refine_extrema(self,
        xlims=None,
        ylims=None,
        xminstep=None,
        yminstep=None,
        minima=True,
        maxima=True,
    ):
        self.debugout("refine extrema: ")
        x = self.x[:]
        y = self.y[:,self.masksensitive]

        if ylims is None:
            ylims = self.ylims
        if ylims is None:
            if self.sameyscale:
                ylims = (asarray(y.ravel().min())[None],asarray(y.ravel().max())[None])
            else:
                ylims = (y.min(axis=0),y.max(axis=0))
        else:
            if not self.sameyscale:
                ylims = (ylims[0][self.masksensitive],ylims[1][self.masksensitive])

        if yminstep is None:
            yminstep = self.yminstep
        if yminstep is None:
            yminstep = (ylims[1] - ylims[0]) * self.precision

        if self.period is not None:
            x = concatenate((x,x[:2]+self.period),axis=0)
            y = concatenate((y,y[:2,:]),axis=0)

        diff = ((y[1:,:] - y[:-1,:])/yminstep).round()
        if minima and maxima:
            def is_extremum(diff0,diff1):
                return (diff0 * diff1 < 0) # | ((diff0 == 0) ^ (diff1 == 0))
        else:
            if minima:
                def is_extremum(diff0,diff1):
                    return ((diff0 < 0) & (diff1 > 0)) # | ((diff0 < 0) & (diff1 >= 0))
            else:
                assert maxima
                def is_extremum(diff0,diff1):
                    return ((diff0 > 0) & (diff1 < 0)) # | ((diff0 > 0) & (diff1 <= 0))

        extrema = zeros(y.shape,bool)
        extrema[1:-1,:] = is_extremum(diff[:-1,:],diff[1:,:])
        if ylims is not None:
            extrema &= (y >= ylims[0])
            extrema &= (y <= ylims[1])
        extrema = extrema.any(1)

        xsplit = (x[1:] + x[:-1]) / 2
        if self.randomize != 0.0:
            xsplit += (rand(len(x)-1) - 0.5) * (x[1:] - x[:-1]) * self.randomize
        self.addpoints(compress(extrema[1:] | extrema[:-1],xsplit),xlims=xlims,xminstep=xminstep)

    def roundoff_extrema(self,xminstep=None):
        self.debugout("roundoff extrema: ")
        xlims = self.xlims
        if xlims is None:
            if self.period is None:
                xlims = (x[0],x[-1])
            else:
                xlims = (0,self.period)
        if xminstep is None:
            xminstep is self.xminstep
        if xminstep is None:
            xminstep = (xlims[1] - xlims[0]) * 1e-3

        x = self.x[self.find_extrema().any(1)]
        self.addpoints(concatenate([x-2*xminstep,x-xminstep,x+xminstep,x+2*xminstep]),xminstep=xminstep)


    def refine_extrema_quad(self,):
        self.debugout("refine extrema quadratically: ")
        x = self.x[:,None]
        y = self.y[:,self.masksensitive]

        if self.period is not None:
            x = concatenate((x,x[:2,:]+self.period),axis=0)
            y = concatenate((y,y[:2,:]),axis=0)

        xdiff = x[1:,:] - x[:-1,:]
        ydiff = y[1:,:] - y[:-1,:]
        slope = ydiff / xdiff
        off = (y[:-1,:]*x[1:,:] - y[1:,:]*x[:-1,:]) / xdiff[:,:]

        xdiff2 = x[2:,:] - x[:-2,:]

        assert allclose(slope * x[:-1,:] + off, y[:-1,:])
        assert allclose(slope * x[1:,:] + off, y[1:,:])

        slopediff = slope[1:,:] - slope[:-1,:]
        curve = slopediff/xdiff2
        xc = (off[:-1,:] - off[1:,:]  + x[:-2,:] * slope[1:,:] - x[2:,:] * slope[:-1,:]) / (2*slopediff)

        # par(x) = curve * (x-xc)**2 + yc
        # line1(x) = (x-x0)*y1/(x1-x0) + (x-x1)*y0/(x0-x1) = x*(y1-y0)/(x1-x0) + (y0*x1 - y1*x0)/(x1-x0)
        #          = slope1 * x + off1
        # line2(x) = (x-x1)*y2/(x2-x1) + (x-x2)*y1/(x1-x2) = x*(y2-y1)/(x2-x1) + (y1*x2 - y2*x1)/(x2-x1)
        #          = slope2 * x + off2
        # par(x) = (x-x0)*line2(x)/(x2-x0) + (x-x2)*line1(x)/(x0-x2)
        #        = ((x-x0)*line2(x) - (x-x2)*line1(x))/(x2-x0)
        #        = (slope2 * (x-x0)*(x+off2/slope2) - slope1*(x-x2)*(x+off1/slope1))/(x2-x0)
        #        = ((slope2-slope1)*x**2 + (off2-x0*slope2 - off1+x2*slope1)*x + yc')/(x2-x0)
        #        = (slope2-slope1)/(x2-x0) * (x + (off2-x0*slope2 - off1+x2*slope1)/(2*(slope2-slope1)))**2 + yc
        #     curve = (slope2 - slope1)/(x2-x0)
        #     xc = (off1 - off2 + x0*slope2 - x2*slope1)/(2*(slope2-slope1))
        #     yc = y1-curve*(x1-xc)**2
        #
        # par(x1) = (x1-x0)*y1/(x2-x0) + (x1-x2)*y1/(x0-x2) = y1*1 ok

        yc = y[1:-1,:] - curve * (x[1:-1,:] - xc)**2
        assert allclose(curve * (x[:-2,:] - xc)**2 + yc, y[:-2,:])
        assert allclose(curve * (x[2:,:] - xc)**2 + yc, y[2:,:])

        extrema = (ydiff[1:,:] * ydiff[:-1,:] < 0) | ((ydiff[1:,:] == 0) ^ (ydiff[:-1,:] == 0))

        assert all((xc - x[:-2,:] > 0)[extrema])
        assert all((x[2:,:] - xc > 0)[extrema])

        self.addpoints(xc[extrema])


    def refine_visible(self,maxangle=pi/100,xlims=None,ylims=None,xminstep=None,yminstep=None):
        self.debugout("refine visible: Nold=%i, "%(len(self.x)))
        x = self.x
        y = self.y[:,self.masksensitive]

        xlims_arg = xlims
        if xlims is None:
            xlims = self.xlims
        if xlims is None:
            if self.period is None:
                xlims = (x[0],x[-1])
            else:
                xlims = (0,self.period)

        if xminstep is None:
            xminstep = self.xminstep
        if xminstep is None:
            xminstep = (xlims[1] - xlims[0]) * self.precision

        if ylims is None:
            ylims = self.ylims
        if ylims is None:
            if self.sameyscale:
                ylims = (asarray(y.ravel().min())[None],asarray(y.ravel().max())[None])
            else:
                ylims = (y.min(axis=0),y.max(axis=0))
        else:
            if not self.sameyscale:
                ylims = (ylims[0][self.masksensitive],ylims[1][self.masksensitive])

        if yminstep is None:
            yminstep = self.yminstep
        if yminstep is None:
            yminstep = (ylims[1] - ylims[0]) * self.precision * 10

        if self.period is not None:
            x = concatenate((
                x[-2:] - self.period,
                x,
                x[:2] + self.period,
            ))
            y = concatenate((
                y[-2:,:],
                y,
                y[:2,:],
            ),axis=0)

        xdiff = x[1:] - x[:-1]
        xdiff2 = x[2:] - x[:-2]
        ydiff = y[1:,:] - y[:-1,:]
        slope = ydiff/xdiff[:,None]
        offset = (y[:-1,:]*x[1:,None] - y[1:,:]*x[:-1,None])/xdiff[:,None]

        y_interpolated = (xdiff[:-1,None]*y[2:,:] + xdiff[1:,None] * y[:-2,:]) / xdiff2[:,None]
        y_interpolated_is_off = (abs(y_interpolated - y[1:-1,:]) > yminstep)

        y_left_extrapolated = x[:-2,None]*slope[1:,:] + offset[1:,:]
        y_right_extrapolated = x[2:,None]*slope[:-1,:] + offset[:-1,:]

        select = zeros(xdiff.shape+y.shape[1:],bool)
        select[:-1,:] |= y_interpolated_is_off
        select[1:,:] |= y_interpolated_is_off
        select[:-1,:] |= (abs(y_left_extrapolated - y[:-2,:]) > yminstep)
        select[1:,:] |= (abs(y_right_extrapolated - y[2:,:]) > yminstep)
        select[:,:] |= isnan(y[:-1,:]) ^ isnan(y[1:,:])

        Nselected = select.any(axis=-1).sum()
        self.debugout("yminstep-selected: %i"%Nselected)

        select &= ((y[1:,:] >= ylims[0][None,:]) | (y[:-1,:] >= ylims[0][None,:]))
        select &= ((y[1:,:] <= ylims[1][None,:]) | (y[:-1,:] <= ylims[1][None,:]))
        select = select.any(1)
        Nylim_cutoff = Nselected - select.sum()
        if Nylim_cutoff>0:
            self.debugout("ylim-cutoff: %i"%Nylim_cutoff)

        select &= xdiff >= 2*xminstep
        Nxminstep_cutoff = Nselected - select.sum() - Nylim_cutoff
        if Nxminstep_cutoff>0:
            self.debugout("xminstep-cutoff0: %i"%Nxminstep_cutoff)

        xsplit = (x[1:] + x[:-1]) / 2
        if self.randomize != 0.0:
            xsplit += (rand(len(x)-1) - 0.5) * (x[1:] - x[:-1]) * self.randomize
        self.addpoints(compress(select,xsplit),xlims=xlims_arg,internalcall=True)


    def reduce_visible(self,xlims=None,ylims=None,xminstep=None,yminstep=None):
        self.debugout("reduce visible: Nold=%i, "%(len(self.x)))
        x = self.x
        y = self.y[:,self.masksensitive]

        xlims_arg = xlims
        if xlims is None:
            xlims = self.xlims
        if xlims is None:
            if self.period is None:
                xlims = (x[0],x[-1])
            else:
                xlims = (0,self.period)

        if xminstep is None:
            xminstep = self.xminstep
        if xminstep is None:
            xminstep = (xlims[1] - xlims[0]) * self.precision

        if ylims is None:
            ylims = self.ylims
        if ylims is None:
            if self.sameyscale:
                ylims = (asarray(y.ravel().min())[None],asarray(y.ravel().max())[None])
            else:
                ylims = (y.min(axis=0),y.max(axis=0))
        else:
            if not self.sameyscale:
                ylims = (ylims[0][self.masksensitive],ylims[1][self.masksensitive])

        if yminstep is None:
            yminstep = self.yminstep
        if yminstep is None:
            yminstep = (ylims[1] - ylims[0]) * self.precision

        if self.period is not None:
            x = concatenate((
                x[-2:] - self.period,
                x,
                x[:2] + self.period,
            ))
            y = concatenate((
                y[-2:,:],
                y,
                y[:2,:],
            ),axis=0)

        xdiff = x[1:] - x[:-1]
        xdiff2 = x[2:] - x[:-2]
        ydiff = y[1:,:] - y[:-1,:]
        slope = ydiff/xdiff[:,None]
        offset = (y[:-1,:]*x[1:,None] - y[1:,:]*x[:-1,None])/xdiff[:,None]

        # linear interpolation: y1 = ((x1-x0)*y2 + (x2-x1)*y0)/(x2-x0)
        y_interpolated = (xdiff[:-1,None]*y[2:,:] + xdiff[1:,None] * y[:-2,:]) / xdiff2[:,None]
        y_interpolated_error = abs(y_interpolated - y[1:-1,:])

        select = ones(self.x.shape,bool)
        select[1:-1] &= (y_interpolated_error > yminstep).any(axis=1)

        self.debugout("interpolatable: %i"%(len(self.x) - select.sum()))

        select[1:(len(select)//2)*2:2] |= ~select[0:(len(select)//2)*2:2]
        select[1:((len(select)-1)//2)*2+1:2] |= ~select[2:((len(select)-1)//2)*2+1:2]

        self.debugout("remove: %i\n"%(len(self.x) - select.sum()))

        self.x = self.x[select]
        self.y = self.y[select,:]


    def refine_valuecut(self,value,xminstep=None):
        shifted = self.y - value
        sgnchange = (shifted[1:,:] * shifted[:-1,:]) < 0
        sel, = sgnchange.any(axis=1).nonzero()
        x0 = self.x[sel][:,None]
        x1 = self.x[sel+1][:,None]
        y0 = self.y[sel,:]
        y1 = self.y[sel+1,:]
        cutx = x0 + (x1-x0)/(y1-y0) * (value - y0)
        selcutx = cutx[sgnchange[sel,:]]
        self.addpoints(selcutx,xminstep=xminstep)
        if xminstep is not None:
            self.addpoints(concatenate((selcutx + xminstep,selcutx - xminstep)),xminstep=xminstep)

    def find_valuecut(self,value,calccutx=True,calcidx=False,calcslope=False,calcslopesign=False):
        x = self.x
        y = self.y
        if self.period is not None:
            x = concatenate((
                x,
                x[:1] + self.period,
            ))
            addyvals = y[:1,:]
            if hasattr(self,'totalpermut'):
                addyvals = addyvals[:,self.totalpermut_inv]
            y = concatenate((
                y,
                addyvals,
            ),axis=0)

        shifted = y - value
        sgn = 2.0 * (shifted >= 0) - 1.0
        sgnleft = sgn[:-1,:]
        sgnright = sgn[1:,:]
        sgnchange = (sgnright - sgnleft)/2.
        sel, = sgnchange.any(axis=1).nonzero()
        x0 = x[sel][:,None]
        x1 = x[sel+1][:,None]
        y0 = y[sel,:]
        y1 = y[sel+1,:]
        slope = (y1-y0)/(x1-x0)
        cutx = x0 + (value - y0) / slope
        res = []
        if calccutx:
            res += [cutx[sgnchange[sel,:]!=0]]
        if calcidx:
            res += [sgnchange[sel,:].nonzero()[1]]
        if calcslope:
            res += [slope[sgnchange[sel,:]!=0]]
        if calcslopesign:
            res += [sgnchange[sel,:][sgnchange[sel,:]!=0]]
        return res


    def interpolate_linear(self,xgrid):
        x = self.x[:]
        y = self.y[:,:]

        if self.period is not None:
            x = concatenate((x[-1:]-self.period,x,x[:1]+self.period),axis=0)
            y = concatenate((y[-1:,:],y,y[:1,:]),axis=0)

        xgrid = asarray(xgrid)

        assert all(xgrid <= x[-1])
        assert all(xgrid >= x[0])

        segm_idx = x.searchsorted(xgrid)
        if segm_idx[0] == 0:
            segm_idx[0] = 1
        segm_idx -= 1
        xdiff = x[1:] - x[:-1]
        ydiff = y[1:] - y[:-1]
        slope = ydiff/xdiff[:,None]
        offset = (y[:-1,:]*x[1:,None] - y[1:,:]*x[:-1,None]) / xdiff[:,None]
        return slope[segm_idx] * xgrid[:,None] + offset[segm_idx]


    def integrate(self):
	res = ((self.x[2:] - self.x[:-2])[:,None] * self.y[1:-1,:]).sum(axis=0)
	if self.period is None:
	    res += (self.x[1] - self.x[0]) * self.y[0,:]
	    res += (self.x[-1] - self.x[-2]) * self.y[-1,:]
	else:
	    res += (self.x[1] - self.x[-1] - self.period) * self.y[0,:]
	    res += (self.x[0] - self.x[-2] - self.period) * self.y[-1,:]
	return 0.5*res


    def average(self):
	if self.period is None:
	    return self.integrate() / (self.x[-1] - self.x[0])
	else:
	    return self.integrate() / (self.period)


if __name__ == '__main__':
    import cnt,chain
    import units

    set_printoptions(precision=2,linewidth=200)

    N = 50
    B = 200

    x = cnt.armchair(N)
    ch = chain.tight_binding_1stNN_graphene(x)
    ch.set_bfield([B * units.Tesla,0,0])

    def do_scan():
        scan = scan_adaptive(
            lambda E: ch.band_energies(E)[2*N-10:2*N+10],
            linspace(2*pi/3*0.9,2*pi/3*1.1,8),
            periodic = False,
#            periodic = True,
        )


        for i in range(10):
            scan.refine_bends(maxquot=2.0)
            print "Now at: %i points"%len(scan.x)
#            scan.refine_extrema()

        globals()['scan'] = scan

    import profile
    profile.run('do_scan()','scanprofile')

    import plot
    for i in range(scan.y.shape[1]):
        plot.plot(scan.x,scan.y[:,i],'-ro')
    plot.show()
