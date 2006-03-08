from calc import *


class scan_adaptive:
    def __init__(self,f,xgrid,periodic=False):
        self.f = f
        self.x = array(sorted(xgrid))
        assert len(self.x.shape) == 1
        assert self.x.shape[0] >= 2
        assert isreal(self.x[0])
        y0 = array(f(self.x[0]))
        self.y = resize(y0,self.x.shape + y0.shape)
        self.y[0,...] = y0
        for i in range(1,len(self.x)):
            self.y[i,...] = array(self.f(self.x[i]))
        self.periodic = periodic
        if periodic:
            assert allclose(self.y[0,...], self.y[-1,...])

    def addpoints(self,newxpoints):
        if self.periodic:
            newxpoints = (newxpoints - self.x[0]) % (self.x[-1] - self.x[0]) + self.x[0]
#	idx,newx = unique1d(concatenate((newxpoints,self.x)).round(),retindx=True)
 #       newy = concatenate((zeros((len(newx),)+self.y.shape[1:]),self.y),axis=0)[idx]

  #      resize(self.y,(len(newx),)+self.y.shape[1:])[idx]

        newxpoints = unique(newxpoints)
        newx = concatenate((self.x,newxpoints))
        newy = resize(self.y,(len(newx),)+self.y.shape[1:])
        for i in range(len(self.x),len(newx)):
            newy[i,...] = array(self.f(newx[i]))
        permut = argsort(newx)
        self.x = take(newx,permut)
        self.y = take(newy,permut,axis=0)

    def divide(self,
        N=2,
        randomize=False, # ToDo
    ):
        newxpoint = []
        for n in range(1,N):
            newxpoints += list(n*self.x[:-1]/N + (N-n)*self.x[:-1]/N)
        self.addpoints(newxpoints)

    def find_extrema(self,):
        diff = self.y[1:,...] - self.y[:-1,...]
        extrema = (diff[1:,...] * diff[:-1,...]) <= 0
        extr_at_boundary = zeros(self.y.shape[1:],bool)
        if self.periodic:
            extr_at_boundary = (diff[0,...] * diff[-1,...]) < 0
        extrema = concatenate(([extr_at_boundary],extrema,[extr_at_boundary]))
        return extrema

    def refine_extrema(self,):
        y = reshape(self.y,(self.y.shape[0],prod(self.y.shape[1:])))
        diff = y[1:,:] - y[:-1,:]
        flatpoints = (diff == 0).any(1)
        extrema = (diff[1:,:] * diff[:-1,:] < 0).any(1)

        extr_at_boundary = False
        if self.periodic:
            extr_at_boundary = (diff[0,:] * diff[-1,:] < 0).any()
        extrema = concatenate(([extr_at_boundary],extrema,[extr_at_boundary]))
        xsplit = (self.x[1:] + self.x[:-1]) / 2
        self.addpoints(compress(extrema[1:] | flatpoints | extrema[:-1],xsplit))

    def refine_extrema_quad(self,):
        x = self.x[:,None]
        y = reshape(self.y,(self.y.shape[0],prod(self.y.shape[1:])))

        if self.periodic:
            x = concatenate((x,x[-1:,:]+(x[1:2,:]-x[0:1,:])),axis=0)
#            x = concatenate((x[0:1,:]-(x[-1:,:]-x[-2:-1,:]),x,x[-1:,:]+(x[1:2,:]-x[0:1,:])),axis=0)
            assert allclose(y[0,:], y[-1,:])
            y = concatenate((y,y[1:2,:]),axis=0)
#            y = concatenate((y[-2:-1,:],y,y[1:2,:]),axis=0)

#        print x

        xdiff = x[1:,:] - x[:-1,:]
#        print xdiff
        xdiff2 = x[2:,:] - x[:-2,:]
        ydiff = y[1:,:] - y[:-1,:]
        slope = ydiff / xdiff
        off = (y[:-1,:]*x[1:,:] - y[1:,:]*x[:-1,:]) / xdiff[:,:]
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
#        print curve * (x[:-2,:] - xc)**2 + yc
#        print y[:-2,:]
        assert allclose(curve * (x[:-2,:] - xc)**2 + yc, y[:-2,:])
        assert allclose(curve * (x[2:,:] - xc)**2 + yc, y[2:,:])

        extrema = (ydiff[1:,:] * ydiff[:-1,:] < 0) | ((ydiff[1:,:] == 0) ^ (ydiff[:-1,:] == 0))

#        print (xc - x[:-2,:])[extrema]
        assert all((xc - x[:-2,:] > 0)[extrema])
#        print (x[2:,:] - xc)[extrema]
        assert all((x[2:,:] - xc > 0)[extrema])

        self.addpoints(xc[extrema])

    def refine_bends(self,maxquot=1.3):
        assert maxquot > 1.0
        y = reshape(self.y,(self.y.shape[0],prod(self.y.shape[1:])))
        diff = y[1:,:] - self.y[:-1,:]
        slope = diff / (self.x[1:,None] - self.x[:-1,None])
        bend = slope[:-1,:] / slope[1:,:]

        issharp = ((bend < 1/maxquot) | (bend > maxquot) | isinf(bend)).any(1)
        issharp_at_boundary = False
        if self.periodic:
            bend_at_boundary = slope[-1,:] / slope[1,:]
            issharp_at_boundary = any((bend_at_boundary < 1/maxquot) | (bend_at_boundary > maxquot) | isinf(bend_at_boundary))
        issharp = concatenate(([issharp_at_boundary],issharp,[issharp_at_boundary]))
        xsplit = (self.x[1:] + self.x[:-1]) / 2

        self.addpoints(compress(issharp[1:] | issharp[:-1],xsplit))

    def refine_visible(self,maxangle=pi/100,xlims=None,ylims=None,xres=1e3,yres=1e2):
        if xlims is None:
            xmin = self.x[0]
            xmax = self.x[-1]
        else:
            xmin,xmax = xlims
        if ylims is None:
            if hasattr(self,'ylims'):
                ymin,ymax = self.ylims
            else:
                ymin = amin(self.y)
                ymax = amax(self.y)
                self.ylims = (ymin,ymax)
        else:
            xmin,xmax = xlims
            self.ylims = ylims

        y = reshape(self.y,(self.y.shape[0],prod(self.y.shape[1:])))

        xdiff = self.x[1:] - self.x[:-1]
        xdiff2 = self.x[2:] - self.x[:-2]

        # linear interpolation: y1 = ((x1-x0)*y2 + (x2-x1)*y0)/(x2-x0)
        y_interpolated = (xdiff[:-1,None]*y[2:,:] + xdiff[1:,None] * y[:-2,:]) / xdiff2[:,None]

        yval_is_off = (abs(y_interpolated - y[1:-1,:]) > (ymax-ymin)/yres).any(1)
        yval_is_off = concatenate(([False],yval_is_off,[False]))
        xdiff_large_enough = xdiff > (xmax-xmin)/xres

        xsplit = (self.x[1:] + self.x[:-1]) / 2

        self.addpoints(compress((yval_is_off[1:] | yval_is_off[:-1]) & xdiff_large_enough,xsplit))

    def refine_sharp_bends_old(self,maxbend=None):
        y = reshape(self.y,(self.y.shape[0],prod(self.y.shape[1:])))
        diff = (y[1:,:] - y[:-1,:]) / (self.x[1:,None] - self.x[:-1,None])
        diff2 = diff[1:,:] - diff[:-1,:]
        maxdiff2 = amax(abs(diff2),axis=-1)
        if maxbend is None:
            maxbend = max(maxdiff2) / 10
        print "maxbend ="
        print maxbend
        sharp_bends = maxdiff2 > maxbend

        sharp_bend_at_boundary = False
        if self.periodic:
            sharp_bend_at_boundary = sometrue(ravel(abs(diff[0,...] - diff[-1,...]) > maxbend))
        sharp_bends = concatenate(([sharp_bend_at_boundary],sharp_bends,[sharp_bend_at_boundary]))
        xsplit = (self.x[1:] + self.x[:-1]) / 2

        self.addpoints(compress(sharp_bends[1:] | sharp_bends[:-1],xsplit))



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




dummy = """
    def refine(self,
        mintol=1e-3,
        division=2,
        minstep_factor=1,
        maxpoints=None,
        randomize=True,  # randomize division points?
    ):
        xminstep = (self.x[-1] - self.x[0]) * mintol * minstep_factor
        if maxpoints is None:
            maxpoints = int(1/mintol)
        if len(self.x) == 2:
            self.divide()








if ylims
    yrange = yrange * ones(1,N);
end

X = 2:(length(x)-1);
xintenv = [ 2*xminstep
            2*xminstep ];
xrelpos = [ 1
            1 ];

yintenv = zeros(2,N);
yexpect = zeros(2,N);
yerr = zeros(2,N);

newx = startgrid(2:(length(startgrid)-1))'; % startgrid was a row-vector'

tol = 256; % start out with a large tol value
while ~isempty(newx)
    newy = zeros(length(newx),N);
    for i=1:length(newx)
        newy(i,:) = feval(fun,newx(i),funparam{:});
    end

    is_new = [zeros(size(x));ones(size(newx))];
    x = [x;newx];
    [x,permut,dummy] = unique(x);

    y = [y;newy]; y=y(permut,:);
    is_new = is_new(permut);
    xintenv = [xintenv;zeros(size(newx))];xintenv = xintenv(permut);
    xrelpos = [xrelpos;zeros(size(newx))];xrelpos = xrelpos(permut);
    yintenv = [yintenv;zeros(size(newy))];yintenv = yintenv(permut,:);
    yexpect = [yexpect;zeros(size(newy))];yexpect = yexpect(permut,:);
    yerr = [yerr;zeros(size(newy))];yerr = yerr(permut,:);

    toadjust = [0;is_new(1:length(x)-2) | is_new(2:length(x)-1) | is_new(3:length(x));0];

    X = find(toadjust);

    xintenv(X) = x(X+1) - x(X-1);
    xrelpos(X) = (x(X) - x(X-1)) ./ xintenv(X);

    for n=1:N
        yintenv(X,n) = y(X+1,n) - y(X-1,n);
        yexpect(X,n) = y(X-1,n) + yintenv(X,n) .* xrelpos(X);
        yerr(X,n) = abs(y(X,n) - yexpect(X,n)) ./ (4 * xrelpos(X) .* (1 - xrelpos(X)));
    end

    newx = [];

    if ~ylims
        for n = 1:N
            [dummy,permut] = sort(-abs(y(:,n)));
%            permut = permut(length(permut):-1:1);
            xsum = 0;
            count = 0;
            while xsum < xrange*0.05;
                count = count + 1;
                xsum = xsum + xintenv(permut(count))/2;
            end
            yrange(n) = abs(y(permut(count)));


%            select = find(y(:,n) ~= 0);
%            if isempty(select)
%                yrange(n) = 0;
%            else
                % a special algorithm that does not react too strongly on
                % sigularities when estimating the y-range of the function
%%                yrange(n) = exp(dot(log(abs(y(select,n))),xintenv(select)/2) / sum(xintenv(select)/2));
%                yrange(n) = dot(abs(y(select,n)),xintenv(select)/2) / sum(xintenv(select)/2);
%            end
        end
    end

    if length(x) >= MAXPOINTS break, end

    while true
        todo = zeros(length(x),1);
        for n=1:N
            todo(X) = todo(X) | (yerr(X,n) > yrange(n) * tol);
        end

        todoleft = todo & (xintenv .* xrelpos > xminstep);
        todoright = todo & (xintenv .* (1-xrelpos) > xminstep);

        newx = [ (x(find(todoleft)-1) + 2*x(find(todoleft)))/3
                 (2*x(find(todoright)) + x(find(todoright)+1))/3 ];

        if ~isempty(newx) break, end

        if tol <= MINTOL break, end

        tol = tol/2;
        if tol < MINTOL
            tol = MINTOL;
        end

        X = 2:(length(x)-1);
    end
    fprintf(1,'done: %i \ttodo: %i \ttol: %f\ttime/sample: %f s\n',length(x),length(newx),tol,etime(clock,start_clock)/length(x));
end

if nargout == 0
    if ~ylims
        for n = 1:N
            ymin(n) = max(min(y(:,n)), -yrange(n)*30);
            ymax(n) = min(max(y(:,n)),  yrange(n)*30);
        end
        ymin = min(ymin);
        ymax = max(ymax);
        ymid = (ymin+ymax)/2;
        yrange = ymax-ymin;
        ymin = ymid - 0.55*yrange;
        ymax = ymid + 0.55*yrange;
    end

    plot(x,y,plotoptions{:})
    set(gca,'XLim',[xmin xmax]);
    set(gca,'YLim',[ymin ymax]);
else
    xres = x;
    yres = y;
end
"""