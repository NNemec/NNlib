from scipy.optimize import leastsq

def fit(x,y,fun,p0,
	return_errors = False):
    if len(p0) == 1:
        def peval(x,p):
	    return fun(x,p)
	def residuals(p,y,x):
    	    return y-peval(x,p)
	assert(x.shape == y.shape)
	assert(fun(x,p0[0]).shape == y.shape)
	xout,cov_x,infodict,ier,mesg = leastsq(residuals,p0,args=(y,x),full_output=True)
	xerr = abs(cov_x[0,0])**.5
	if return_errors:
    	    return (xout,),(xerr,)
	else:
	    return (xout,)
	
    else:
        def peval(x,p):
	    return fun(x,*p)
	def residuals(p,y,x):
    	    return y-peval(x,p)
	assert(x.shape == y.shape)
	assert(fun(x,*p0).shape == y.shape)
	xout,cov_x,infodict,ier,mesg = leastsq(residuals,p0,args=(y,x),full_output=True)
	xerr = []
        for i in range(len(xout)):
    	    xerr += [abs(cov_x[i,i])**.5]
	if return_errors:
    	    return xout,xerr
	else:
	    return xout

if __name__ == '__main__':
    from scipy import *

    def myfun(x,a,b,c,d):
        return a*exp(-b*x**2-c)+b*x+c

    (a,b,c,d) = (2.,3.,-2.,1.)
    x = linspace(0,10,10)
    y = myfun(x,a,b,c,d)

    (a,b,c,d) = (0.,0.,0.,0.)
    (a,b,c,d) = fit(x,y,myfun,(a,b,c,d))
    print a,b,c,d

    import pylab
    pylab.plot(x,y)
    pylab.plot(x,myfun(x,a,b,c,d))
    pylab.show()
