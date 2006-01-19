#!/usr/bin/env python

from calc import *
from units import *
from param import param
import xyz


def SDOS_armchair_analytic(
    E, # in units of gamma
    N, # chirality of tube: (N,N)
):
    Gs = 0.j
    for j in range(1,2*N+1):
        X = E**2-sin(j*pi/N)**2
        if X < 0:
            continue
        cos_q = -.5*(cos(j*pi/N) + sign(E)*(X)**.5)
        if cos_q**2 > 1:
            continue
        print cos_q
        sin_q = (1-cos_q**2)**.5
        Gs += E*(.5+1j*sin_q/(2*cos_q+cos(j*pi/N)))
    return -Gs.imag/pi


def armchair(N):
    CC_distance = param.GRAPHENE_CC_DISTANCE

    period = 3**.5 * CC_distance
    r = N*2*1.5*CC_distance / (2*pi)

    res = xyz.chain((0,0,period))
    res.radius = r
    for n in range(N):
        res.atoms.extend([
	    xyz.atom('C',(r*cos(2*pi*(6*n  )/(6*N)),r*sin(2*pi*(6*n  )/(6*N)),period/2)),
	    xyz.atom('C',(r*cos(2*pi*(6*n+1)/(6*N)),r*sin(2*pi*(6*n+1)/(6*N)),0)),
	    xyz.atom('C',(r*cos(2*pi*(6*n+3)/(6*N)),r*sin(2*pi*(6*n+3)/(6*N)),0)),
	    xyz.atom('C',(r*cos(2*pi*(6*n+4)/(6*N)),r*sin(2*pi*(6*n+4)/(6*N)),period/2)),
	])

    return res


def zigzag(N):
    CC_distance = param.GRAPHENE_CC_DISTANCE

    period = 3*CC_distance
    r = N * 3.**.5 * CC_distance / (2*pi)

    res = xyz.chain((0,0,period))
    res.radius = r
    for n in range(N):
        res.atoms.extend([
    	    xyz.atom('C',(r*cos(2*pi*(2*n  )/(2*N)),r*sin(2*pi*(2*n  )/(2*N)),period/6)),
    	    xyz.atom('C',(r*cos(2*pi*(2*n  )/(2*N)),r*sin(2*pi*(2*n  )/(2*N)),period/2)),
    	    xyz.atom('C',(r*cos(2*pi*(2*n+1)/(2*N)),r*sin(2*pi*(2*n+1)/(2*N)),0)),
    	    xyz.atom('C',(r*cos(2*pi*(2*n+1)/(2*N)),r*sin(2*pi*(2*n+1)/(2*N)),period/1.5)),
	])

    return res

def gcd(a,b):
    if a>b:
        a,b = b,a
    while a != 0:
        a,b = b%a,a
    return b

def is_metallic(M,N):
    assert M >= 0
    assert N >= 0
    if M==0:
        M,N=N,0
        assert M > 0
    return ((M-N)%3 == 0)

A_plaquette = a**2 * 3**.5 * 1.5

def radius(M,N):
    assert M >= 0
    assert N >= 0
    if M==0:
        M,N=N,0
        assert M > 0
    CC_distance = param.GRAPHENE_CC_DISTANCE
    metallic = ((M-N)%3 == 0)
    circumference = CC_distance * sqrt(3.0) * sqrt(M**2 + N**2 + M*N)
    return circumference/(2*pi)

def A_section(M,N):
    return radius(M,N)**2 * pi

def Natoms(M,N):
    assert M >= 0
    assert N >= 0
    assert M+N > 0
    multiple_perp = gcd((M+2*N),(2*M+N))
    Nplaquettes = 2 * (M**2 + N**2 + M*N)/multiple_perp
    return Nplaquettes * 2

def chiral(M,N):
    assert M >= 0
    assert N >= 0
    if M==0:
        M,N=N,0
        assert M > 0
    CC_distance = param.GRAPHENE_CC_DISTANCE
    r = radius(M,N)
    multiple = gcd(M,N)
    multiple_perp = gcd((M+2*N),(2*M+N))
    M_perp = (M+2*N) / multiple_perp
    N_perp = -(2*M+N) / multiple_perp
    period = CC_distance * sqrt(3.0) * sqrt(M_perp**2 + N_perp**2 + M_perp*N_perp);
    Nplaquettes = 2 * (M**2 + N**2 + M*N)/multiple_perp
    Natoms = Nplaquettes * 2

    Nlines = (M+N)/multiple
    start = [0]*Nlines
    stop = [0]*Nlines
    for l in range(Nlines):
        start[l] = -(l*N/(M+N))
    for l in range(Nlines):
        stop[l] = - N_perp + start[(l-(M_perp+N_perp))%Nlines] - N/multiple*((l-(M_perp+N_perp))/Nlines)

#    assert (sum(stop) - sum(start)) * multiple == Nplaquettes

    """
    X*a = R*(Ma+Nb)+L*((M+2N)a-(2M+N)b)

    b: 0 = RN-2LM-LN
       R = 2LM/N+L
    a: X = (2LM/N+L)M+LM+2LN

    L:=N
    R = 2M+N
    X = 2MM+2NM+2NN
    """

    dphi_a = 2*pi * (2*M+N)/(2*(M*M + M*N + N*N))
    dphi_b = 2*pi * (M+2*N)/(2*(M*M + M*N + N*N))
    dphi_c = dphi_a - dphi_b

    dz_a = CC_distance * sqrt(3) * sqrt((2*M+N)**2 + (M+2*N)**2 - (2*M+N)*(M+2*N)) * N / (2*(M*M + M*N + N*N))
    dz_b = - CC_distance * sqrt(3) * sqrt((2*M+N)**2 + (M+2*N)**2 - (2*M+N)*(M+2*N)) * M / (2*(M*M + M*N + N*N))
    dz_c = dz_a - dz_b

#    for i in sorted(locals().keys()):
#        print i+":\t"+str(locals()[i])

    res = xyz.chain((0,0,period))
    res.radius = r

    n = 0
    for l in range(Nlines):
        for p in range(start[l],stop[l]):
            for m in range(multiple):
                phi_A = dphi_a * l + dphi_c * p + 2*pi*m/multiple
                z_A = dz_a * l + dz_c * p
                phi_B = phi_A + (dphi_a+dphi_b)/3
                z_B = z_A + (dz_a+dz_b)/3
                res.atoms.extend([
		    xyz.atom('C',(cos(phi_A)*r,sin(phi_A)*r,z_A)),
		    xyz.atom('C',(cos(phi_B)*r,sin(phi_B)*r,z_B)),
		])
    assert len(res.atoms) == Natoms

    return res


def swcnt(V):
    assert len(V)==2
    if V[0] == V[1]:
        return armchair(V[0])
    elif V[1] == 0:
        return zigzag(V[0])
    elif V[0] == 0:
        return zigzag(V[1])
    else:
	return chiral(V[0],V[1])


if __name__ == "__main__":
    if False:
        from plot import *
        N = 6

        E = linspace(-8.5,8.5,200)*eV

        SDOS = array([
            SDOS_armchair_analytic(e/(2.66*eV),N)/(2.66*eV*N)
            for e in E
        ])

        def integral(x,y):
            return .5*sum((x[1:]-x[:-1])*(y[1:]+y[:-1]))

        print "integral: %g"%integral(E,SDOS)

        plot(E,SDOS)
        show()
    else:
#        param.setdefaults()
        zigzag(10).write_xyz_file('cnt-test-zigzag.xyz')
        armchair(10).write_xyz_file('cnt-test-armchair.xyz')
        chiral(5,4).write_xyz_file('cnt-test-chiral.xyz')

