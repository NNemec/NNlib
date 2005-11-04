from calc import *

import bfield as bf
import xyz

from param import param

class chain:
    def __init__(self,xyz_chain,H_int_B0,H_hop_B0):
        assert type(H_int_B0) is type(Matrix(()))
        assert type(H_hop_B0) is type(Matrix(()))
        assert isinstance(xyz_chain,xyz.chain)
        self.xyz = xyz_chain
        self.xyz_shifted = xyz_chain.shift(xyz_chain.period)
        self.N_atoms = len(self.xyz.atoms)
        assert(shape(H_int_B0) == (self.N_atoms,self.N_atoms))
        assert(shape(H_hop_B0) == (self.N_atoms,self.N_atoms))
        self.H_int_B0 = H_int_B0
        self.H_hop_B0 = H_hop_B0
        self.bfield = None
        self.set_bfield(0)
        self.energy = None

    def set_bfield(self,bfield):
        if bfield != self.bfield:
            self.bfield = bfield
            if self.bfield == 0:
                self.H_int = self.H_int_B0
                self.H_hop = self.H_hop_B0
            else:
                self.H_int = bf.calc_H_int(self.bfield,self.H_int_B0,self.xyz)
                self.H_hop = bf.calc_H_hop(self.bfield,self.H_hop_B0,self.xyz,self.xyz_shifted)
            self.cache = {}

    def set_energy(self,energy):
        if energy is None:
            assert self.energy is not None
        else:
            self.energy = energy
        if self.energy in self.cache:
            (self._G_bulk,self._Gs_L,self._Gs_R) = self.cache[self.energy]
        else:
            self._do_calc_lopez_sancho()
            self.cache[self.energy] = (self._G_bulk,self._Gs_L,self._Gs_R)

    def _do_calc_lopez_sancho(self):
        # ToDo: find documentation (Lopez-Sancho)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)*Matrix(eye(size(self.H_int,0)))

        # alpha = energy*S_hop - chain.H_hop;
        # beta = energy*chain.S_hop' - chain.H_hop';
        alpha = - self.H_hop
        beta = - adj(self.H_hop)
        epsilon = E - self.H_int
        epsilon_L = epsilon
        epsilon_R = epsilon

        i = 0;
        EPS = param.LOPEZ_SANCHO_EPSILON;
        while norm(alpha,1) > EPS or norm(beta,1) > EPS:
            gamma = inv(epsilon)
            temp_agb = alpha*gamma*beta
            temp_bga = beta*gamma*alpha
            alpha = alpha*gamma*alpha
            beta = beta*gamma*beta
            epsilon = epsilon - temp_agb - temp_bga
            epsilon_L = epsilon_L - temp_bga
            epsilon_R = epsilon_R - temp_agb

            i = i+1;
            if i > 100:
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

    def band_energies(self,k):
        H_eff = self.H_int + exp(1j*k)*self.H_hop + exp(-1j*k)*adj(self.H_hop)
        return array(sorted(list(real(eigvals(H_eff)))))

    def DOS(self,energy=None):
        return -1./pi*imag(trace(self.G_bulk(energy)))/self.N_atoms

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

def tight_binding_1stNN_graphene(xyz_chain):
    (H_int,H_hop,N) = _tight_binding_1stNN_graphene_H(xyz_chain)
    return chain(xyz_chain,H_int,H_hop)

def tight_binding_dwcnt_triozon(xyz_tube_A,xyz_tube_B):
    # based on the parametrization described in
    # doi:10.1103/PhysRevB.64.121401
    
    def interwall_hopping(pos_a,pos_b):
        d = norm(pos_b-pos_a);
        cos_theta = vdot(pos_a[:2],pos_b[:2])/(norm(pos_a[:2])*norm(pos_b[:2]));
        return -param.TRIOZON_BETA * cos_theta * exp((param.TRIOZON_A - d)/param.TRIOZON_DELTA);

    xyz_chain = xyz.merge(xyz_tube_A,xyz_tube_B)

    (H_int_A,H_hop_A,N_A) = _tight_binding_1stNN_graphene_H(xyz_tube_A)
    (H_int_B,H_hop_B,N_B) = _tight_binding_1stNN_graphene_H(xyz_tube_B)
    H_int = Matrix(zeros((N_A+N_B,N_A+N_B),'D'))
    H_int[:N_A,:N_A] = H_int_A
    H_int[N_A:,N_A:] = H_int_B
    H_hop = Matrix(zeros((N_A+N_B,N_A+N_B),'D'))
    H_hop[:N_A,:N_A] = H_hop_A
    H_hop[N_A:,N_A:] = H_hop_B

    for b in range(N_B):
        pos_b = xyz_tube_B.atoms[b].pos
        pos_b_hop = pos_b + xyz_chain.period
        for a in range(N_A):
            pos_a = xyz_tube_A.atoms[a].pos
            pos_a_hop = pos_a + xyz_chain.period
            H_int[a,N_A+b] = interwall_hopping(pos_a,pos_b)
            H_int[N_A+b,a] = H_int[a,N_A+b]
            H_hop[a,N_A+b] = interwall_hopping(pos_a,pos_b_hop)
            H_hop[N_A+b,a] = interwall_hopping(pos_b,pos_a_hop)

    return chain(xyz_chain,H_int,H_hop)


if __name__ == "__main__":
    x = xyz.square_ladder(3)
    x3 = x.multiply(3)
    ch = chain(x3)

    print ch.band_energies(pi/2)
