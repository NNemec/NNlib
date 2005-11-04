from calc import *

class scan_adaptive:
    def __init__(self,f,xgrid,periodic=False):
        self.f = f
        self.x = array(sorted(xgrid))
        assert len(self.x.shape) == 1
        assert self.x.shape[0] >= 2
#        assert self.x.typecode in ('d','D','f','F')
        y0 = array(f(self.x[0]))
#        tc = y0.typecode
#        assert tc in ('d','D','f','F')
        self.y = resize(y0,self.x.shape + y0.shape)
#        self.y[0,...] = y0
        for i in range(1,len(self.x)):
            self.y[i,...] = array(self.f(self.x[i]))
        self.periodic = periodic
        if periodic:
            assert self.y[0,...] == self.y[-1,...]

    def addpoints(self,newxpoints):
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

    def findextrema(self,):
        diff = self.y[1:,...] - self.y[:-1,...]
        flatpoints = sometrue(reshape(
            diff == 0,
            (diff.shape[0],prod(diff.shape[1:])),
        ),axis=-1)
        extrema = sometrue(reshape(
            (diff[1:,...] * diff[:-1,...]) < 0,
            (diff.shape[0]-1,prod(diff.shape[1:])),
        ),axis=-1)

        extr_at_boundary = False
        if self.periodic:
            extr_at_boundary = sometrue(ravel((diff[0,...] * diff[-1,...]) < 0))
        extrema = concatenate(([extr_at_boundary],extrema > 0,[extr_at_boundary]))
        xsplit = (self.x[1:] + self.x[:-1]) / 2
        self.addpoints(compress(extrema[1:] | flatpoints | extrema[:-1],xsplit))




if __name__ == '__main__':
    import xyz,chain
    from param import param
    param.setdefaults()
    x = xyz.armchair(3)
    ch = chain.tight_binding_1stNN_graphene(x)

    def do_scan():
        scan = scan_adaptive(
            ch.band_energies,
            linspace(0,2*pi,50),
            periodic = True,
        )


        for i in range(10):
            scan.findextrema()

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