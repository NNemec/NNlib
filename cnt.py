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

    period = c_[0,0,sqrt(3)*CC_distance]
    circumference = N*2*1.5*CC_distance
    radius = circumference/(2*pi)
    angle = 2*pi/(6*N)

    rot = Matrix(eye(3))
    rot[0,0] = cos(angle)
    rot[1,1] = cos(angle)
    rot[0,1] = sin(angle)
    rot[1,0] = -sin(angle)
    rot2 = rot*rot
    rot3 = rot2*rot
    rot4 = rot2*rot2
    rot6 = rot3*rot3

    at = xyz.atom('C',c_[radius,0,0])

    cell = xyz.structure()
    cell.atoms = [
        at.shift(period/2),
        at.rotate(rot),
        at.rotate(rot3),
        at.rotate(rot4).shift(period/2)
    ]

    res = xyz.chain(period)
    for i in range(N):
        res.atoms.extend(cell.atoms)
        cell = cell.rotate(rot6)
    res.radius = radius

    return res

def zigzag(N):
    CC_distance = param.GRAPHENE_CC_DISTANCE

    period = c_[0,0,3*CC_distance]
    circumference = N*sqrt(3)*CC_distance
    radius = circumference/(2*pi)
    angle = pi/N

    rot = Matrix(eye(3))
    rot[0,0] = cos(angle)
    rot[1,1] = cos(angle)
    rot[0,1] = sin(angle)
    rot[1,0] = -sin(angle)
    rot2 = rot*rot

    at = xyz.atom('C',c_[radius,0,0])

    cell = xyz.structure()
    cell.atoms = [
        at.shift(period/6),
        at.shift(period/2),
        at.rotate(rot),
        at.rotate(rot).shift(period/1.5)
    ]

    res = xyz.chain(period)
    for i in range(N):
        res.atoms.extend(cell.atoms)
        cell = cell.rotate(rot2)
    res.radius = radius

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
    Natoms = Nplaquettes * 2

def chiral(M,N):
    assert M >= 0
    assert N >= 0
    if M==0:
        M,N=N,0
        assert M > 0
    CC_distance = param.GRAPHENE_CC_DISTANCE
    radius = radius(M,N)
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

    atoms = [None] * Natoms

    n = 0
    for l in range(Nlines):
        for p in range(start[l],stop[l]):
            for m in range(multiple):
                phi_A = dphi_a * l + dphi_c * p + 2*pi*m/multiple
                z_A = dz_a * l + dz_c * p
                phi_B = phi_A + (dphi_a+dphi_b)/3
                z_B = z_A + (dz_a+dz_b)/3
                atoms[n] = xyz.atom('C',c_[cos(phi_A)*radius,sin(phi_A)*radius,z_A])
                atoms[n+1] = xyz.atom('C',c_[cos(phi_B)*radius,sin(phi_B)*radius,z_B])
                n += 2

    res = xyz.chain(period)
    for i in range(Natoms):
        res.atoms = atoms
    res.radius = radius

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
        param.setdefaults()
        ch = chiral(5,4)
        ch.write_xyz_file('cnt-test.xyz')

