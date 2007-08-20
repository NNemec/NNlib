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
        self.N_orbitals = H_B0[0].shape[0]
        for h_b0 in H_B0:
            assert type(h_b0) is matrix
            assert h_b0.shape == (self.N_orbitals,self.N_orbitals)
        self.H = H_B0
        self.H_B0 = H_B0

        if xyz_chain is not None:
            assert isinstance(xyz_chain,xyz.chain)
            self.xyz = xyz_chain
            self.xyz_shifted = [ xyz_chain.shift(xyz_chain.period * i) for i in range(1,len(H_B0)) ]
            self.bfield = array((0,0,0))

        if S is not None:
            for s in S:
                assert type(s) is matrix
                assert s.shape == (self.N_orbitals,self.N_orbitals)
            self.S = S

        self.energy = None
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
        assert not hasattr(self,'S')
        assert len(self.H) == 2

        # ToDo: find documentation (Lopez-Sancho)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA) * Matrix(eye(self.N_orbitals))

        # alpha = energy*S_hop - chain.H_hop;
        # beta = energy*chain.S_hop' - chain.H_hop';
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
        def adjsum(a):
            return a + a.H
        for i in range(1,len(self.H)):
            res += adjsum(exp(1j*k*i)*self.H[i])
        return res

    def S_eff(self,k):
        if not hasattr(self,'S'):
            return None
        res = self.S[0].copy()
        def adjsum(a):
            return a + a.H
        for i in range(1,len(self.S)):
            res += adjsum(exp(1j*k*i)*self.S[i])
        return res

    def band_energies(self,k):
        if hasattr(self,'S'):
            return array(sorted(list(real(scipy.linalg.eigvals(self.H_eff(k),self.S_eff(k))))))
        else:
            return array(sorted(list(real(eigvalsh(self.H_eff(k))))))

    def DOS(self,energy=None):
        return -1./pi*imag(trace(self.G_bulk(energy)))/self.N_orbitals

    def LDOS(self,energy=None):
        return -1./pi*imag(diag(self.G_bulk(energy)))

    def SDOS_L(self,energy=None):
        return -1./pi*imag(trace(self.Gs_L(energy)))/self.N_orbitals

    def SDOS_R(self,energy=None):
        return -1./pi*imag(trace(self.Gs_L(energy)))/self.N_orbitals

    def transmission(self,energy=None):
        self.set_energy(energy)
        E = (self.energy+1j*param.LOPEZ_SANCHO_ETA)*Matrix(eye(self.N_orbitals))
        Sigma_L = adj(self.H[1])*self.Gs_L()*self.H[1]
        Sigma_R = self.H[1]*self.Gs_R()*adj(self.H[1])
        Gamma_L = 1j*(Sigma_L-adj(Sigma_L))
        Gamma_R = 1j*(Sigma_R-adj(Sigma_R))
        Gc = inv(E-self.H[0]-Sigma_L-Sigma_R)
        return real(trace(Gamma_L*Gc*Gamma_R*adj(Gc)))

    def multiply(self,N):
        xyz = None
        A = self.N_orbitals
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

if __name__ == "__main__":
    import cnt
    import tightbinding
    x = cnt.armchair(20)
    ch = tightbinding.tight_binding_1stNN_graphene(x)

    a = ch.Gs_L(energy=0.5)
    b = ch.Gs_L(energy=1.0)
