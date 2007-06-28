#!/usr/bin/env python

class square_interp:
    def __init__(self,xy0,xy1,xy2):
        # yl(x) = (y1*(x-x0) + y0*(x1-x))/(x1-x0)
        #       = (y1-y0)/(x1-x0) * x + (x1*y0-x0*y1)/(x1-x0)
        #       = al * x + bl
        # yr(x) = (y2*(x-x1) + y1*(x2-x))/(x2-x1)
        #       = (y2-y1)/(x2-x1) * x + (x2*y1-x1*y2)/(x2-x1)
        #       = ar * x + br
        #
        # y(x) = (yr(x)*(x-x0) + yl(x)*(x2-x))/(x2-x0)
        #      = ((ar * x + br)*(x-x0) + (al * x + bl)*(x2-x))/(x2-x0)
        #      = ((ar-al) * x^2 + ((br-bl)+(al*x2-ar*x0)) * x + (bl*x2 - br*x0))/(x2-x0)
        #      = a*x^2 + b*x + c
        #      = a*(x + b/2a)^2 - b^2/4a + c
        #      = a(x-x_min)^2 + y_min

        x0,y0 = xy0
        x1,y1 = xy1
        x2,y2 = xy2
        al = (y1-y0) * 1./(x1-x0)
        ar = (y2-y1) * 1./(x2-x1)
        bl = (x1*y0-x0*y1) * 1./(x1-x0)
        br = (x2*y1-x1*y2) * 1./(x2-x1)
        b = ((br-bl)+(al*x2-ar*x0)) * 1./(x2-x0)
        a = (ar-al) * 1./(x2-x0)
        c = (bl*x2 - br*x0) * 1./(x2-x0)

        x_min = - b / (2*a)
        y_min = c - b**2 / (4*a)

        self.__dict__.update(locals())

    def yl(self,x):
        return self.al * x + self.bl

    def yr(self,x):
        return self.ar * x + self.br

    def y(self,x):
        return self.a * x**2 + self.b * x + self.c

if __name__ == '__main__':
    sq = square_interp((0,1),(1,0),(2,2))
    print sq.yl(0)
    print sq.yl(1)
    print sq.yr(1)
    print sq.yr(2)
    print sq.y(0)
    print sq.y(1)
    print sq.y(2)
#    print sq.__dict__
