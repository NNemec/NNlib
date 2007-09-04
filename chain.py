from calc import *
import xyz
from param import param
from units import *

param.createdefault("LOPEZ_SANCHO_ETA", 1e-5*eV)
param.createdefault("LOPEZ_SANCHO_EPSILON", 1e-4*eV)
param.createdefault("LOPEZ_SANCHO_MAXSTEPS", 100)

class chain:
    def __init__(self,H_B0,xyz_chain=None,do_cache=True,S=None):
        assert type(H_B0) is list
        N = H_B0[0].shape[0]
        self.N_orbitals = N
        for h_b0 in H_B0:
            assert type(h_b0) is matrix
            assert h_b0.shape == (N,N)
	    assert h_b0.dtype == H_B0[0].dtype
        self.H = H_B0
        self.H_B0 = H_B0

        if xyz_chain is not None:
            assert isinstance(xyz_chain,xyz.chain)
            self.xyz = xyz_chain
            self.xyz_shifted = [ xyz_chain.shift(xyz_chain.period * i) for i in range(1,len(H_B0)) ]
            self.bfield = array((0,0,0))

        self.nonorthogonal = (S is not None)

        if self.nonorthogonal:
            for s in S:
                assert type(s) is matrix
                assert s.shape == (N,N)
	        assert s.dtype == S[0].dtype
            self.S = S
        else:
            self.S = [ matrix(eye(N)) ] + [ matrix(zeros((N,N))) ] * (len(H_B0)-1)

        self.energy = None
        self.Sigma = 0.0
        if do_cache:
            self.cache = {}
        (self._G_bulk,self._Gs_L,self._Gs_R) = (None,None,None)

    def set_bfield(self,bfield):
        import bfield as bf
        assert hasattr(self,'xyz')
        assert len(self.xyz.atoms) == self.N_orbitals
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

    def set_energy(self,energy,Sigma=0.0):
        if energy is None:
            assert self.energy is not None
        else:
            if energy != self.energy or any(Sigma != self.Sigma):
                (self._G_bulk,self._Gs_L,self._Gs_R) = (None,None,None)
            self.energy = energy
	self.Sigma = Sigma
        if self._G_bulk is None:
            if hasattr(self,'cache') and self.energy in self.cache and all(Sigma == 0.0):
                (self._G_bulk,self._Gs_L,self._Gs_R) = self.cache[self.energy]
            else:
                self._do_calc_lopez_sancho(Sigma=Sigma)
                if hasattr(self,'cache') and all(Sigma == 0.0):
                    self.cache[self.energy] = (self._G_bulk,self._Gs_L,self._Gs_R)

    def _do_calc_lopez_sancho(self,Sigma=0.0):
        assert len(self.H) == 2

        # ToDo: find documentation (Lopez-Sancho)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)

        epsilon = E*self.S[0] - self.H[0] - Sigma

        if self.nonorthogonal:
            alpha = E*self.S[1] - self.H[1]
            beta = E*self.S[1].H - self.H[1].H
        else:
            alpha = - self.H[1]
            beta = - self.H[1].H

        epsilon_L = epsilon
        epsilon_R = epsilon

        i = 0;
        EPS = param.LOPEZ_SANCHO_EPSILON;
        while abs(alpha).A.sum() + abs(beta).A.sum() > EPS:
            gamma = epsilon.I
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

        G_bulk = epsilon.I
        Gs_L = epsilon_L.I
        Gs_R = epsilon_R.I

        self._G_bulk = G_bulk
        self._Gs_L = Gs_L
        self._Gs_R = Gs_R

    def Gs_L(self,energy=None,Sigma=0.0):
        self.set_energy(energy,Sigma)
        return self._Gs_L

    def Gs_R(self,energy=None,Sigma=0.0):
        self.set_energy(energy,Sigma)
        return self._Gs_R

    def G_bulk(self,energy=None,Sigma=0.0):
        self.set_energy(energy,Sigma)
        return self._G_bulk

    def H_eff(self,k):
        res = self.H[0] + 0.0j
        def adjsum(a):
            return a + a.H
        for i in range(1,len(self.H)):
            res += adjsum(exp(1j*k*i)*self.H[i])
        return res

    def S_eff(self,k):
        res = self.S[0] + 0.0j
        if self.nonorthogonal:
            def adjsum(a):
                return a + a.H
            for i in range(1,len(self.S)):
                res += adjsum(exp(1j*k*i)*self.S[i])
        return res

    def band_energies(self,k):
        if self.nonorthogonal:
#            X = self.S_eff(k).I * self.H_eff(k)
#            return array(sorted(list(real(scipy.linalg.eigvals(X)))))
            return array(sorted(list(real(scipy.linalg.eigvals(self.H_eff(k),self.S_eff(k))))))
        else:
            return array(sorted(list(real(eigvalsh(self.H_eff(k))))))

    def DOS(self,energy=None,Sigma=0.0):
        return -1./pi*imag(trace(self.G_bulk(energy,Sigma)))/self.N_orbitals

    def LDOS(self,energy=None,Sigma=0.0):
        return -1./pi*imag(diag(self.G_bulk(energy,Sigma)))

    def SDOS_L(self,energy=None,Sigma=0.0):
        return -1./pi*imag(trace(self.Gs_L(energy,Sigma)))/self.N_orbitals

    def SDOS_R(self,energy=None,Sigma=0.0):
        return -1./pi*imag(trace(self.Gs_L(energy,Sigma)))/self.N_orbitals

    def transmission(self,energy=None, Sigma=0.0):
        assert len(self.H) == 2
        self.set_energy(energy, Sigma=Sigma)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)
        if self.nonorthogonal:
            ES_H_1 = E*self.S[1] - self.H[1]
            ESh_Hh_1 = E*self.S[1].H - self.H[1].H
        else:
            ES_H_1 = - self.H[1]
            EaS_aH_1 = - self.H[1].H
        Sigma_L = ESh_Hh_1 * self.Gs_L() * ES_H_1
        Sigma_R = ES_H_1 * self.Gs_R() * ESh_Hh_1
        Gamma_L = 1j*(Sigma_L - Sigma_L.H)
        Gamma_R = 1j*(Sigma_R - Sigma_R.H)
        Gc = (E*self.S[0]-self.H[0]-Sigma-Sigma_L-Sigma_R).I
        return real(trace(Gamma_L * Gc * Gamma_R * Gc.H))

    def multiply(self,N = None):
        if N == None:
            N = len(self.H_B0) - 1
        xyz = None
        A = self.N_orbitals
        if hasattr(self,'xyz'):
            xyz = self.xyz.multiply(N)
        assert len(self.H_B0) <= N+1
        H = [ matrix(zeros((N*A,N*A),self.H_B0[0].dtype)) for i in range(2) ]
        for n in range(N):
            H[0][n*A:(n+1)*A,n*A:(n+1)*A] = self.H_B0[0]
        for i in range(1,len(self.H_B0)):
            for n in range(i):
                H[1][(n-i+N)*A:(n-i+N+1)*A,n*A:(n+1)*A] = self.H_B0[i]
            for n in range(i,N):
                H[0][(n-i)*A:(n-i+1)*A,n*A:(n+1)*A] = self.H_B0[i]
                H[0][n*A:(n+1)*A,(n-i)*A:(n-i+1)*A] = self.H_B0[i].H
        if self.nonorthogonal:
            S = [ matrix(zeros((N*A,N*A),self.S[0].dtype)) for i in range(2) ]
            for n in range(N):
                S[0][n*A:(n+1)*A,n*A:(n+1)*A] = self.S[0]
            for i in range(1,len(self.S)):
                for n in range(i):
                    S[1][(n-i+N)*A:(n-i+N+1)*A,n*A:(n+1)*A] = self.S[i]
                for n in range(i,N):
                    S[0][(n-i)*A:(n-i+1)*A,n*A:(n+1)*A] = self.S[i]
                    S[0][n*A:(n+1)*A,(n-i)*A:(n-i+1)*A] = self.S[i].H
        else:
            S = None
        return chain(H,xyz,S=S)

def square_ladder(N,gamma,gamma_perp=None,do_cache=True):
    if gamma_perp == None:
        gamma_perp = gamma

    H = [ matrix(zeros((N,N))) for i in range(2) ]

    for n in range(1,N):
        H[0][n-1,n] = -gamma_perp
        H[0][n,n-1] = -gamma_perp

    for n in range(N):
        H[1][n,n] = -gamma

    return chain(H,do_cache=do_cache)

def linchain(gamma,do_cache=True):
    return square_ladder(N=1,gamma=gamma,do_cache=do_cache)

if __name__ == "__main__":
    import cnt
    import tightbinding
    x = cnt.armchair(20)
    ch = tightbinding.tight_binding_1stNN_graphene(x)

    a = ch.Gs_L(energy=0.5)
    b = ch.Gs_L(energy=1.0)
