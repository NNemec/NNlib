from calc import *

import xyz

from param import param
from units import *

class sheet:
    def __init__(self,H_B0,xyz_sheet=None,S=None):
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

        if S is not None:
            assert type(S) is dict
            assert S.keys() == H_B0.keys()
            for k in S:
                assert type(S[k]) is matrix
                assert S[k].shape == (self.N_atoms,self.N_atoms)
            self.S = S

        if xyz_sheet is not None:
            assert isinstance(xyz_sheet,xyz.sheet)
            assert len(xyz_sheet.atoms) == self.N_atoms
            self.xyz = xyz_sheet
            self.latticecoords = None
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
                if self.latticecoords is None:
                    self.latticecoords = [ (dot(at.pos,self.xyz.rzp[0]) , dot(at.pos,self.xyz.rzp[1]) ) for at in self.xyz.atoms ]

                self.H = deepcopy(self.H_B0)
                assert all(bfield[:2] == 0.) # bfield perpendicular to sheet
                assert self.xyz.period[0][2] == 0.
                assert self.xyz.period[1][2] == 0.

                def phase(s,d,shift):
                    # linear gauge:
                    # phase = (d[1] - s[1]) * (d[0] + s[0]) * .5

                    shiftx,shifty = shift

                    si = int(s[0] // 1)
                    sr = s[0] % 1
                    sx = si+sr
                    sy = s[1]
                    di = int(d[0] // 1) + shiftx
                    dr = d[0] % 1
                    dx = di+dr
                    dy = d[1] + shifty

                    # forced gauge:
                    if dx == sx:
                        phase = (dy - sy) * (dr)

                    else:
                        Aavg = ( (di - si) * .5 + (dr - sr) * (dr + sr) * .5 ) / (dx - sx)
                        phase = (dy - sy) * Aavg
                        if di > si:
                            for i in range(si,di):
                                x = i+1
                                y = x * (dy-sy)/(dx-sx) + (dx*sy-sx*dy)/(dx-sx)
                                phase -= y
                        elif di < si:
                            for i in range(di,si):
                                x = i+1
                                y = x * (dy-sy)/(dx-sx) + (dx*sy-sx*dy)/(dx-sx)
                                phase += y

                    return phase%1


                if self.phase is None:
                    self.phase = {}

                    for i0,i1 in self.H:
                        self.phase[i0,i1] = zeros(self.H[i0,i1].shape)

                        for n,m in array(asarray(self.H[i0,i1]).nonzero()).transpose():
                            self.phase[i0,i1][n,m] = phase(self.latticecoords[n],self.latticecoords[m],(i0,i1))

                for i0,i1 in self.H:
                    for n,m in array(asarray(self.H[i0,i1]).nonzero()).transpose():
                        self.H[i0,i1][n,m] *= exp(1j*2*pi*self.phase[i0,i1][n,m]*Nflux)


    def H_eff(self,ka):
        res = self.H[0,0] + 0.0
        def adjsum(a):
            return a + a.H
        for i0,i1 in self.H:
            if (i0,i1) != (0,0):
                res += adjsum(exp(1j*dot(ka,[i0,i1]))*self.H[i0,i1])
        return res

    def S_eff(self,ka):
        if not hasattr(self,'S'):
            return None
        res = self.S[0,0].copy()
        def adjsum(a):
            return a + a.H
        for i0,i1 in self.S:
            if (i0,i1) != (0,0):
                res += adjsum(exp(1j*dot(ka,[i0,i1]))*self.S[i0,i1])
        return res

    def band_energies(self,ka):
        if hasattr(self,'S'):
            return array(sorted(list(real(scipy.linalg.eigvals(self.H_eff(ka),self.S_eff(ka))))))
        else:
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
