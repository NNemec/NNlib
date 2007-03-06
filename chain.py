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
    def __init__(self,H_B0,xyz_chain=None,do_cache=True,S=None):
        assert type(H_B0) is list
        self.N_atoms = H_B0[0].shape[0]
        for h_b0 in H_B0:
            assert type(h_b0) is matrix
            assert h_b0.shape == (self.N_atoms,self.N_atoms)
        self.H = H_B0
        self.H_B0 = H_B0

        if xyz_chain is not None and len(xyz_chain.atoms) == self.N_atoms:
            assert isinstance(xyz_chain,xyz.chain)
            assert len(xyz_chain.atoms) == self.N_atoms
            self.xyz = xyz_chain
            self.xyz_shifted = [ xyz_chain.shift(xyz_chain.period * i) for i in range(1,len(H_B0)) ]
            self.bfield = array((0,0,0))

        if S is not None:
            for s in S:
                assert type(s) is matrix
                assert s.shape == (self.N_atoms,self.N_atoms)
            self.S = S

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
                self.H = self.H_B0
            else:
                self.H = (
                    [ bf.calc_H_int(self.bfield,self.H_B0[0],self.xyz) ]
                    +
                    [ bf.calc_H_hop(self.bfield,self.H_B0[i],self.xyz,self.xyz_shifted[i-1]) for i in range(1,len(self.H_B0)) ]
                )
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
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA) * Matrix(eye(self.N_atoms))

        # alpha = energy*S_hop - chain.H_hop;
        # beta = energy*chain.S_hop' - chain.H_hop';
        assert len(self.H) == 2
        assert not hasattr(self,'S')
        alpha = - self.H[1]
        beta = - adj(self.H[1])
        epsilon = E - self.H[0]
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

        G_bulk = inv(epsilon)
        Gs_L = inv(epsilon_L)
        Gs_R = inv(epsilon_R)

        self._G_bulk = G_bulk
        self._Gs_L = Gs_L
	self._Gs_R = Gs_R

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
        res = self.H[0] + 0.0
        for i in range(1,len(self.H)):
            res += exp(1j*k*i)*self.H[i] + exp(-1j*k*i)*adj(self.H[i])
        return res

    def S_eff(self,k):
        if not hasattr(self,'S'):
            return None
        res = self.S[0] + 0.0
        for i in range(1,len(self.S)):
            res += exp(1j*k*i)*self.S[i] + exp(-1j*k*i)*adj(self.S[i])
        return res

    def band_energies(self,k):
        if hasattr(self,'S'):
            return array(sorted(list(real(scipy.linalg.eigvals(self.H_eff(k),self.S_eff(k))))))
        else:
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
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)*Matrix(eye(self.N_atoms))
        Sigma_L = adj(self.H[1])*self.Gs_L()*self.H[1]
        Sigma_R = self.H[1]*self.Gs_R()*adj(self.H[1])
        Gamma_L = 1j*(Sigma_L-adj(Sigma_L))
        Gamma_R = 1j*(Sigma_R-adj(Sigma_R))
        Gc = inv(E-self.H[0]-Sigma_L-Sigma_R)
        return real(trace(Gamma_L*Gc*Gamma_R*adj(Gc)))

    def multiply(self,N):
        xyz = None
        A = self.N_atoms
        if hasattr(self,'xyz'):
            xyz = self.xyz.multiply(N)
        assert len(self.H_B0) <= N+1
        H = [ Matrix(zeros((N*A,N*A),'D')) for i in range(2) ]
        for n in range(N):
            H[0][n*A:(n+1)*A,n*A:(n+1)*A] = self.H_B0[0]
        for i in range(1,len(self.H_B0)):
            for n in range(i):
                H[1][(n-i+N)*A:(n-i+N+1)*A,n*A:(n+1)*A] = self.H_B0[i]
            for n in range(i,N):
                H[0][(n-i)*A:(n-i+1)*A,n*A:(n+1)*A] = self.H_B0[i]
                H[0][n*A:(n+1)*A,(n-i)*A:(n-i+1)*A] = adj(self.H_B0[i])
        if hasattr(self,'S'):
            S = [ Matrix(zeros((N*A,N*A),'D')) for i in range(2) ]
            for n in range(N):
                S[0][n*A:(n+1)*A,n*A:(n+1)*A] = self.S[0]
            for i in range(1,len(self.S)):
                for n in range(i):
                    S[1][(n-i+N)*A:(n-i+N+1)*A,n*A:(n+1)*A] = self.S[i]
                for n in range(i,N):
                    S[0][(n-i)*A:(n-i+1)*A,n*A:(n+1)*A] = self.S[i]
                    S[0][n*A:(n+1)*A,(n-i)*A:(n-i+1)*A] = adj(self.S[i])
        else:
            S = None
        return chain(H,xyz,S=S)

def square_ladder(N,gamma,gamma_perp=None,do_cache=True):
    if gamma_perp == None:
	gamma_perp = gamma
	
    H = [ Matrix(zeros((N,N),'D')) for i in range(2) ]

    for n in range(1,N):
        H[0][n-1,n] = -gamma_perp
        H[0][n,n-1] = -gamma_perp

    for n in range(N):
        H[1][n,n] = -gamma

    return chain(H,do_cache=do_cache)

def linchain(gamma,do_cache=True):
    return square_ladder(N=1,gamma=gamma,do_cache=do_cache)

def tight_binding_graphene_1stNN(xyz_chain,do_cache=True):
    N = len(xyz_chain.atoms)
    H = [ Matrix(zeros((N,N),'D')) for i in range(2) ]
    maxdist = param.GRAPHENE_CC_DISTANCE * 1.1
    gamma = param.GRAPHENE_1STNN_HOPPING

    for i in range(N):
        for j in range(i+1,N):
            if norm(xyz_chain.atoms[i].pos - xyz_chain.atoms[j].pos) < maxdist:
                H[0][i,j] = -gamma
                H[0][j,i] = -gamma
    for i in range(N):
        for j in range(N):
            if norm(xyz_chain.atoms[i].pos - (xyz_chain.atoms[j].pos+xyz_chain.period)) < maxdist:
                H[1][i,j] = -gamma

    return chain(H,xyz_chain,do_cache=do_cache)

def tight_binding_triozon(xyz,do_cache=True,graphite=False):
    # based on the parametrization described in
    # doi:10.1103/PhysRevB.64.121401

    CC_DIST = param.GRAPHENE_CC_DISTANCE
    NN_HOP = param.GRAPHENE_1STNN_HOPPING

    TRIO_CUTOFF = param.TRIOZON_CUTOFF
    Z_CUTOFF = param.TRIOZON_Z_CUTOFF
    BETA = param.TRIOZON_BETA
    A = param.TRIOZON_A
    DELTA = param.TRIOZON_DELTA

    if graphite:
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
    else: # MWCNT
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

    at = xyz.atoms
    period = xyz.period

    Natoms = len(at)
    H = [ Matrix(zeros((Natoms,Natoms),'D')) ]

    for i in range(Natoms):
        for j in range(i+1,Natoms):
            hop = hopping(at[i].pos,at[j].pos)
            if hop != 0.0:
                H[0][i,j] = hop
                H[0][j,i] = conj(hop)

    for n in range(1,100):
        xyz_shifted = xyz.shift(period*n)
        at_sh = xyz_shifted.atoms

        h = Matrix(zeros((Natoms,Natoms),'D'))
        nonzero = False
        for i in range(Natoms):
            for j in range(Natoms):
                hop = hopping(at[i].pos,at_sh[j].pos)
                if hop != 0.0:
                    nonzero = True
                    h[i,j] = hop
        if nonzero:
            H.append(h)
        else:
            break

    assert n < 99

    return chain(H,xyz,do_cache=do_cache)

if __name__ == "__main__":
    import cnt
    x = cnt.armchair(20)
    ch = tight_binding_1stNN_graphene(x)

    a = ch.Gs_L(energy=0.5)
    b = ch.Gs_L(energy=1.0)
