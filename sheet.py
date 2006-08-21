from calc import *

import xyz
import cnt
import chain

from param import param
from units import *

class sheet:
    def __init__(self,H_B0,xyz_sheet=None):
        assert type(H_B0) is dict
        assert (0,0) in H_B0
        assert type(H_B0[0,0]) is matrix
        self.N_atoms = H_B0[0,0].shape[0]
        self.imax = 0
        for i0,i1 in H_B0:
            assert type(i0) is int
            assert type(i1) is int
            assert i0 >= 0
            assert i0 > 0 or i1 >= 0
            self.imax = max(self.imax,i0,i1,-i1)
            assert type(H_B0[i0,i1]) is matrix
            assert H_B0[i0,i1].shape == (self.N_atoms,self.N_atoms)

        self.H_B0 = H_B0
        self.H = H_B0

        if xyz_sheet is not None:
            assert isinstance(xyz_sheet,xyz.sheet)
            assert len(xyz_sheet.atoms) == self.N_atoms
            self.xyz = xyz_sheet
            self.xyz_shifted = None
            self.bfield = array((0,0,0))
            self.area = cross(self.xyz.period[0],self.xyz.period[1])[2]
            self.phase = None



    def set_bfield(self,bfield):
        import bfield as bf
        from copy import deepcopy

        assert hasattr(self,'xyz')

        bfield = asarray(bfield)
        flux = self.area*bfield[2]
        Nflux = round(flux/Phi0)
        bfield[2] = Nflux * Phi0 / self.area

        if any(bfield != self.bfield):
            self.bfield = bfield

            if sum(array(bfield)**2) == 0:
                self.H = self.H_B0
            else:
                if self.xyz_shifted is None:
                    self.xyz_shifted = {}
                    for i0,i1 in self.H_B0:
                        self.xyz_shifted[i0,i1] = self.xyz.shift(i0 * self.xyz.period[0] + i1 * self.xyz.period[1])

                self.H = deepcopy(self.H_B0)
                assert all(bfield[:2] == 0.) # bfield perpendicular to sheet
                assert self.xyz.period[0][2] == 0.
                assert self.xyz.period[1][2] == 0.

                def phase(s,d,):#Nflux):
                    # linear gauge:
                    # phase = (d[1] - s[1]) * (d[0] + s[0]) * .5
#                   print s,d,
                    sx,sy = s
                    dx,dy = d

                    # forced gauge:
                    if dx == sx:
                        phase = (dy - sy) * (dx%1)

                    else:
                        Aavg = ( (dx // 1 - sx // 1) * .5 + (dx % 1 - sx % 1) * (dx % 1 + sx % 1) * .5 ) / (dx - sx)
                        phase = (dy - sy) * Aavg
                        if dx > sx:
                            for i in range(int(sx//1),int(dx//1)):
                                x = i+1
                                y = x * (dy-sy)/(dx-sx) + (dx*sy-sx*dy)/(dx-sx)
                                phase -= y%1
                        else:
                            for i in range(int(dx//1),int(sx//1)):
                                x = i+1
                                y = x * (dy-sy)/(dx-sx) + (dx*sy-sx*dy)/(dx-sx)
                                phase += y%1

#                   print phase

                    return phase

#                    return exp(1j*2*pi*phase*Nflux)


                if self.phase is None:
                    self.phase = {}

                    for i0,i1 in self.H:
                        self.phase[i0,i1] = zeros(self.H[i0,i1].shape)

                        for n,m in array(asarray(self.H[i0,i1]).nonzero()).transpose():
                            pos = self.xyz.atoms[n].pos
                            rs = [ dot(pos,self.xyz.rzp[0]) , dot(pos,self.xyz.rzp[1]) ]
                            pos_shifted = self.xyz_shifted[i0,i1].atoms[m].pos
                            rd = [ dot(pos_shifted,self.xyz.rzp[0]) , dot(pos_shifted,self.xyz.rzp[1]) ]

                            self.phase[i0,i1][n,m] = phase(rs,rd)

                for i0,i1 in self.H:
#                   H = asarray(self.H[i0,i1])
#                   H[H.nonzero()] *= exp(1j*2*pi*self.phase[i0,i1][H.nonzero()]*Nflux)
#                   self.H[i0,i1] = asmatrix(H)

#                   self.H[i0,i1] = asmatrix(asarray(self.H[i0,i1]) * exp(1j*2*pi*self.phase[i0,i1]*Nflux))

                    for n,m in array(asarray(self.H[i0,i1]).nonzero()).transpose():
                        self.H[i0,i1][n,m] *= exp(1j*2*pi*self.phase[i0,i1][n,m]*Nflux)

#                        pos = self.xyz.atoms[n].pos
#                        rs = [ dot(pos,self.xyz.rzp[0]) , dot(pos,self.xyz.rzp[1]) ]
#                        pos_shifted = self.xyz_shifted[i0,i1].atoms[m].pos
#                        rd = [ dot(pos_shifted,self.xyz.rzp[0]) , dot(pos_shifted,self.xyz.rzp[1]) ]
#
#                        self.H[i0,i1][n,m] *= phase(rs,rd,Nflux)
#                   print

    def H_eff(self,ka):
        res = self.H[0,0] + 0.0
        def adjsum(a):
            return a + adj(a)
        for i0,i1 in self.H:
            if (i0,i1) != (0,0):
                res += adjsum(exp(1j*dot(ka,[i0,i1]))*self.H[i0,i1])
        return res

    def band_energies(self,ka):
        return array(sorted(list(real(eigvalsh(self.H_eff(ka))))))

    def multiply(self,N0,N1):
        xyz = None
        Nat = self.N_atoms
        if hasattr(self,'xyz'):
            xyz = self.xyz.multiply(N0,N1)
        H = {}
        H[0,0] = Matrix(zeros((N0*N1*Nat,N0*N1*Nat),'D'))
        for n0 in range(N0):
            for n1 in range(N1):
                H[0,0][
                    (n0*N1+n1)*Nat:(n0*N1+n1+1)*Nat,
                    (n0*N1+n1)*Nat:(n0*N1+n1+1)*Nat,
                ] = self.H[0,0]
                for i0,i1 in self.H:
                    if (i0,i1) == (0,0):
                        continue

                    i0n = (n0+i0) // N0
                    i1n = (n1+i1) // N1
                    if i0n > 0 or (i0n == 0 and i1n >= 0):
                        if (i0n,i1n) not in H:
                            H[i0n,i1n] = Matrix(zeros((N0*N1*Nat,N0*N1*Nat),'D'))

                        H[i0n,i1n][
                            (n0*N1+n1)*Nat:(n0*N1+n1+1)*Nat,
                            ((n0+i0)%N0*N1+(n1+i1)%N1)*Nat:((n0+i0)%N0*N1+(n1+i1)%N1+1)*Nat,
                        ] = self.H[i0,i1]

                    i0n = (n0-i0) // N0
                    i1n = (n1-i1) // N1
                    if i0n > 0 or (i0n == 0 and i1n >= 0):
                        if (i0n,i1n) not in H:
                            H[i0n,i1n] = Matrix(zeros((N0*N1*Nat,N0*N1*Nat),'D'))

                        H[i0n,i1n][
                            (n0*N1+n1)*Nat:(n0*N1+n1+1)*Nat,
                            ((n0-i0)%N0*N1+(n1-i1)%N1)*Nat:((n0-i0)%N0*N1+(n1-i1)%N1+1)*Nat,
                        ] = adj(self.H[i0,i1])

        return sheet(H,xyz)

def square_lattice(gamma):
    H = {}
    H[0,0] = Matrix([[0j]])
    H[0,1] = Matrix([[-gamma + 0j]])
    H[1,0] = Matrix([[-gamma + 0j]])
    return sheet(H)

def graphene(gamma):
    H = {}
    H[0,0] = Matrix([[0j,-gamma],[-gamma,0j]])
    H[0,1] = Matrix([[0j,-gamma],[0j,0j]])
    H[1,0] = Matrix([[0j,0j],[0j,-gamma]])
    return sheet(H)

def tight_binding_1stNN_graphene(xyz_sheet):
    N = len(xyz_sheet.atoms)
    H = {}
    H[0,0] = Matrix(zeros((N,N),'D'))

    maxdist = param.GRAPHENE_CC_DISTANCE * 1.1
    gamma = param.GRAPHENE_1STNN_HOPPING

    for i in range(N):
        for j in range(i+1,N):
            if norm(xyz_sheet.atoms[i].pos - xyz_sheet.atoms[j].pos) < maxdist:
                H[0,0][i,j] = -gamma
                H[0,0][j,i] = -gamma

    for i0,i1 in [(0,1),(1,1),(1,0),(1,-1)]:
        shift = i0 * xyz_sheet.period[0] + i1 * xyz_sheet.period[1]
 #       print "shift: ",shift
        H_hop = Matrix(zeros((N,N),'D'))
        nonzero = False

        for i in range(N):
            for j in range(N):
#                print xyz_sheet.atoms[i].pos - (xyz_sheet.atoms[j].pos + shift)
                if norm(xyz_sheet.atoms[i].pos - (xyz_sheet.atoms[j].pos + shift)) < maxdist:
                    H_hop[i,j] = -gamma
                    nonzero = True
        if nonzero:
            H[i0,i1] = H_hop
 #       print

    return sheet(H,xyz_sheet)


def tight_binding_graphite_triozon(xyz_sheet_A,xyz_sheet_B):
    # based on the parametrization described in
    # doi:10.1103/PhysRevB.64.121401

    CC_DIST = param.GRAPHENE_CC_DISTANCE
    NN_HOP = param.GRAPHENE_1STNN_HOPPING

    TRIO_CUTOFF = param.TRIOZON_CUTOFF
    Z_CUTOFF = param.TRIOZON_Z_CUTOFF
    BETA = param.TRIOZON_BETA
    A = param.TRIOZON_A
    DELTA = param.TRIOZON_DELTA

    def hopping(pos_a,pos_b):
#        if abs(pos_a[2] - pos_b[2]) > Z_CUTOFF:
#            return 0.0
        if abs(pos_a[2] - pos_b[2]) < CC_DIST*0.1:
            if norm(pos_a - pos_b) < CC_DIST*1.1:
                return -NN_HOP
        else:
            d = norm(pos_b-pos_a);
            if d < TRIO_CUTOFF:
                return -BETA * exp((A - d)/DELTA);
        return 0.0

    x = xyz.merge(xyz_sheet_A,xyz_sheet_B)
    at = x.atoms
    period = x.period

    Natoms = len(x.atoms)
    H = {}
    H[0,0] = Matrix(zeros((Natoms,Natoms),'D'))

    for i in range(Natoms):
        for j in range(i+1,Natoms):
            hop = hopping(at[i].pos,at[j].pos)
            if hop != 0.0:
                H[0,0][i,j] = hop
                H[0,0][j,i] = conj(hop)

    for i0 in range(5):
        for i1 in range(-5,5):
            shift = i0 * period[0] + i1 * period[1]
            if norm(shift) > Z_CUTOFF + norm(period[0]) + norm(period[1]):
                continue
            if i0 == 0 and i1 <= 0:
                continue
            h_hop = Matrix(zeros((Natoms,Natoms),'D'))
            nonzero = False

            x_shifted = x.shift(shift)
            at_sh = x_shifted.atoms

            for i in range(Natoms):
                for j in range(Natoms):
                    hop = hopping(at[i].pos,at_sh[j].pos)
                    if hop != 0.0:
                        h_hop[i,j] = hop
                        nonzero = True
            if nonzero:
                H[i0,i1] = h_hop

    return sheet(H,x)
