from calc import *

import xyz
import cnt

from param import param
from units import *

param.createdefault("LOPEZ_SANCHO_ETA", 1e-5*eV)
param.createdefault("LOPEZ_SANCHO_EPSILON", 1e-4*eV)
param.createdefault("LOPEZ_SANCHO_MAXSTEPS", 100)

param.createdefault("TRIOZON_BETA", param.GRAPHENE_1STNN_HOPPING / 8)
param.createdefault("TRIOZON_A", 0.334 * nm)
param.createdefault("TRIOZON_DELTA", 0.045 * nm)
param.createdefault("TRIOZON_CUTOFF", param.TRIOZON_A+5*param.TRIOZON_DELTA)
param.createdefault("TRIOZON_Z_CUTOFF", (param.TRIOZON_CUTOFF**2 - (0.95*param.GRAPHITE_INTERLAYER_DISTANCE)**2)**0.5)

class chain:
    def __init__(self,H_int_B0,H_hop_B0,xyz_chain=None,do_cache=True):
        assert type(H_int_B0) is matrix
        self.N_atoms = H_int_B0.shape[0]
        assert H_int_B0.shape == (self.N_atoms,self.N_atoms)
        if type(H_hop_B0) is matrix:
            H_hop_B0 = [H_hop_B0]
        for hhop in H_hop_B0:
            assert type(hhop) is matrix
            assert hhop.shape == (self.N_atoms,self.N_atoms)
        self.H_int = H_int_B0
        self.H_hop = H_hop_B0
        self.H_int_B0 = H_int_B0
        self.H_hop_B0 = H_hop_B0

        if xyz_chain is not None:
            assert isinstance(xyz_chain,xyz.chain)
            assert len(xyz_chain.atoms) == self.N_atoms
            self.xyz = xyz_chain
            self.xyz_shifted = [ xyz_chain.shift(xyz_chain.period * (i+1)) for i in range(len(H_hop_B0)) ]
            self.bfield = array((0,0,0))

        self.energy = None
        if do_cache:
            self.cache = {}
        (self._G_bulk,self._Gs_L,self._Gs_R) = (None,None,None)


    def set_bfield(self,bfield):
        import bfield as bf
        assert hasattr(self,'xyz')
        if any(bfield != self.bfield):
            self.bfield = bfield
            if sum(array(self.bfield)**2) == 0:
                self.H_int = self.H_int_B0
                self.H_hop = self.H_hop_B0
            else:
                self.H_int = bf.calc_H_int(self.bfield,self.H_int_B0,self.xyz)
                self.H_hop = [ bf.calc_H_hop(self.bfield,self.H_hop_B0[i],self.xyz,self.xyz_shifted[i]) for i in range(len(self.H_hop_B0)) ]
            if hasattr(self,'cache'):
                self.cache = {}
            (self._G_bulk,self._Gs_L,self._Gs_R) = (None,None,None)

    def set_energy(self,energy):
        if energy is None:
            assert self.energy is not None
        else:
            if energy != self.energy:
                (self._G_bulk,self._Gs_L,self._Gs_R) = (None,None,None)
            self.energy = energy
        if self._G_bulk is None:
            if hasattr(self,'cache') and self.energy in self.cache:
                (self._G_bulk,self._Gs_L,self._Gs_R) = self.cache[self.energy]
            else:
                self._do_calc_lopez_sancho()
                if hasattr(self,'cache'):
                    self.cache[self.energy] = (self._G_bulk,self._Gs_L,self._Gs_R)

    def _do_calc_lopez_sancho(self):
        # ToDo: find documentation (Lopez-Sancho)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)*Matrix(eye(size(self.H_int,0)))

        # alpha = energy*S_hop - chain.H_hop;
        # beta = energy*chain.S_hop' - chain.H_hop';
        assert len(self.H_hop) == 1
        alpha = - self.H_hop[0]
        beta = - adj(self.H_hop[0])
        epsilon = E - self.H_int
        epsilon_L = epsilon
        epsilon_R = epsilon

        i = 0;
        EPS = param.LOPEZ_SANCHO_EPSILON;
        while abs(alpha).A.sum() + abs(beta).A.sum() > EPS:
            gamma = inv(epsilon)
            temp_agb = alpha*gamma*beta
            temp_bga = beta*gamma*alpha
            alpha = alpha*gamma*alpha
            beta = beta*gamma*beta
            epsilon = epsilon - temp_agb - temp_bga
            epsilon_L = epsilon_L - temp_bga
            epsilon_R = epsilon_R - temp_agb

            i = i+1;
            if i > param.LOPEZ_SANCHO_MAXSTEPS:
                raise "Lopez Sancho does not converge"

        self._G_bulk = inv(epsilon)
        self._Gs_L = inv(epsilon_L)
        self._Gs_R = inv(epsilon_R)

    def Gs_L(self,energy=None):
        self.set_energy(energy)
        return self._Gs_L

    def Gs_R(self,energy=None):
        self.set_energy(energy)
        return self._Gs_R

    def G_bulk(self,energy=None):
        self.set_energy(energy)
        return self._G_bulk

    def H_eff(self,k):
        res = self.H_int + 0.0
        for i in range(len(self.H_hop)):
            res += exp(1j*k*(i+1))*self.H_hop[i] + exp(-1j*k*(i+1))*adj(self.H_hop[i])
        return res

    def band_energies(self,k):
        return array(sorted(list(real(eigvalsh(self.H_eff(k))))))

    def DOS(self,energy=None):
        return -1./pi*imag(trace(self.G_bulk(energy)))/self.N_atoms

    def LDOS(self,energy=None):
        return -1./pi*imag(diag(self.G_bulk(energy)))

    def SDOS_L(self,energy=None):
        return -1./pi*imag(trace(self.Gs_L(energy)))/self.N_atoms

    def SDOS_R(self,energy=None):
        return -1./pi*imag(trace(self.Gs_L(energy)))/self.N_atoms

    def transmission(self,energy=None):
        self.set_energy(energy)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)*Matrix(eye(size(self.H_int,0)))
        Sigma_L = adj(self.H_hop)*self.Gs_L()*self.H_hop
        Sigma_R = self.H_hop*self.Gs_R()*adj(self.H_hop)
        Gamma_L = 1j*(Sigma_L-adj(Sigma_L))
        Gamma_R = 1j*(Sigma_R-adj(Sigma_R))
        Gc = inv(E-self.H_int-Sigma_L-Sigma_R)
        return real(trace(Gamma_L*Gc*Gamma_R*adj(Gc)))

    def multiply(self,N):
        xyz = None
        A = self.N_atoms
        if hasattr(self,'xyz'):
            xyz = self.xyz.multiply(N)
        H_int = Matrix(zeros((N*A,N*A),'D'))
        assert len(H_hop_B0) == 1
        for n in range(N):
            H_int[n*A:(n+1)*A,n*A:(n+1)*A] = self.H_int_B0
        for n in range(1,N):
            H_int[(n-1)*A:n*A,n*A:(n+1)*A] = self.H_hop_B0[0]
            H_int[n*A:(n+1)*A,(n-1)*A:n*A] = adj(self.H_hop_B0[0])
        H_hop = [Matrix(zeros((N*A,N*A),'D'))]
        H_hop[0][(N-1)*A:N*A,0:A] = self.H_hop_B0[0]
        return chain(H_int,H_hop,xyz)

def square_ladder(N,gamma,do_cache=True):
    H_int = Matrix(zeros((N,N),'D'))
    H_hop = Matrix(zeros((N,N),'D'))

    for n in range(1,N):
        H_int[n-1,n] = -gamma
        H_int[n,n-1] = -gamma

    for n in range(N):
        H_hop[n,n] = -gamma

    return chain(H_int,H_hop,do_cache)

def linchain(gamma,do_cache=True):
    return square_ladder(N=1,gamma=gamma,do_cache=do_cache)

def _tight_binding_1stNN_graphene_H(xyz_chain):
    N = len(xyz_chain.atoms)
    H_int = Matrix(zeros((N,N),'D'))
    H_hop = Matrix(zeros((N,N),'D'))
    maxdist = param.GRAPHENE_CC_DISTANCE * 1.1
    gamma = param.GRAPHENE_1STNN_HOPPING

    for i in range(N):
        for j in range(i+1,N):
            if norm(xyz_chain.atoms[i].pos - xyz_chain.atoms[j].pos) < maxdist:
                H_int[i,j] = -gamma
                H_int[j,i] = -gamma
    for i in range(N):
        for j in range(N):
            if norm(xyz_chain.atoms[i].pos - (xyz_chain.atoms[j].pos+xyz_chain.period)) < maxdist:
                H_hop[i,j] = -gamma
    return (H_int,H_hop,N)

def tight_binding_1stNN_graphene(xyz_chain,do_cache=True):
    (H_int,H_hop,N) = _tight_binding_1stNN_graphene_H(xyz_chain)
    return chain(H_int,H_hop,xyz_chain,do_cache=do_cache)

def tight_binding_dwcnt_triozon(xyz_tube_A,xyz_tube_B,do_cache=True):
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
        if abs(pos_a[2] - pos_b[2]) > Z_CUTOFF:
            return 0.0
        elif abs(norm(pos_a[:2])-norm(pos_b[:2])) < CC_DIST*0.1:
            if norm(pos_a - pos_b) < CC_DIST*1.1:
                return -NN_HOP
        else:
            d = norm(pos_b-pos_a);
            if d < TRIO_CUTOFF:
                cos_theta = vdot(pos_a[:2],pos_b[:2])/(norm(pos_a[:2])*norm(pos_b[:2]));
                return -BETA * cos_theta * exp((A - d)/DELTA);
        return 0.0

    x = xyz.merge(xyz_tube_A,xyz_tube_B)
    at = x.atoms
    period = x.period

    Natoms = len(x.atoms)
    H_int = Matrix(zeros((Natoms,Natoms),'D'))

    for i in range(Natoms):
        for j in range(i+1,Natoms):
            hop = hopping(at[i].pos,at[j].pos)
            if hop != 0.0:
                H_int[i,j] = hop
                H_int[j,i] = conj(hop)

    H_hop = []
    for n in range(20):
#    for n in range(1+int(Z_CUTOFF/period[2])):
#    for n in [0]:
        x_shifted = x.shift(period*(n+1))
        at_sh = x_shifted.atoms

        h_hop = Matrix(zeros((Natoms,Natoms),'D'))
        nonzero = False
        for i in range(Natoms):
            for j in range(Natoms):
                hop = hopping(at[i].pos,at_sh[j].pos)
                if hop != 0.0:
                    nonzero = True
                    h_hop[i,j] = hop
        if nonzero:
            H_hop.append(h_hop)
        else:
            break

    return chain(H_int,H_hop,x,do_cache=do_cache)

def tight_binding_graphite_triozon(xyz_tube_A,xyz_tube_B,do_cache=True):
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
        if abs(pos_a[2] - pos_b[2]) > Z_CUTOFF:
            return 0.0
        elif abs(pos_a[1]-pos_b[1]) < CC_DIST*0.1:
            if norm(pos_a - pos_b) < CC_DIST*1.1:
                return -NN_HOP
        else:
            d = norm(pos_b-pos_a);
            if d < TRIO_CUTOFF:
                return -BETA * exp((A - d)/DELTA);
        return 0.0

    x = xyz.merge(xyz_tube_A,xyz_tube_B)
    at = x.atoms
    period = x.period

    Natoms = len(x.atoms)
    H_int = Matrix(zeros((Natoms,Natoms),'D'))

    for i in range(Natoms):
        for j in range(i+1,Natoms):
            hop = hopping(at[i].pos,at[j].pos)
            if hop != 0.0:
                H_int[i,j] = hop
                H_int[j,i] = conj(hop)

    H_hop = []
    for n in range(20):
#    for n in range(1+int(Z_CUTOFF/period[2])):
#    for n in [0]:
        x_shifted = x.shift(period*(n+1))
        at_sh = x_shifted.atoms

        h_hop = Matrix(zeros((Natoms,Natoms),'D'))
        nonzero = False
        for i in range(Natoms):
            for j in range(Natoms):
                hop = hopping(at[i].pos,at_sh[j].pos)
                if hop != 0.0:
                    nonzero = True
                    h_hop[i,j] = hop
        if nonzero:
            H_hop.append(h_hop)
        else:
            break

    return chain(H_int,H_hop,x,do_cache=do_cache)


if __name__ == "__main__":
    import cnt
    x = cnt.armchair(20)
    ch = tight_binding_1stNN_graphene(x)

    a = ch.Gs_L(energy=0.5)
    b = ch.Gs_L(energy=1.0)
