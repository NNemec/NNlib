from calc import *

from copy import deepcopy
from param import param

param.createdefault("DISORDER_TYPE", "diagonal")
# param.DISORDER_TYPE = "hopping"

param.createdefault("DO_CALC_CHANNELS", False)
param.createdefault("NO_DISORDER_IN_CONTACTS", False)

class conductor:
    def __init__(self,xyz,H_int,H_hop):
        N = [ H.shape[0] for H in H_int ]
#        N = [ len(x.atoms) for x in xyz ]
        assert len(H_int)==len(N)
        assert len(H_hop)==len(N)-1
        for i in range(len(N)):
            assert type(H_int[i]) == type(Matrix(()))
            assert shape(H_int[i]) == (N[i],N[i])
        for i in range(len(N)-1):
            assert type(H_hop[i]) == type(Matrix(()))
            assert shape(H_hop[i]) == (N[i],N[i+1])
        self.N = N
	if xyz is not None:
	    assert len(xyz) == len(N)
	    for i in range(len(N)):
		assert len(xyz[i].atoms) == N[i]
            self.xyz = xyz
        self.H_int_B0 = H_int
        self.H_hop_B0 = H_hop
        self.bfield = array((0,0,0))
        self.H_int = H_int
        self.H_hop = H_hop
        self.disorder = 0
        self.disorder_seed = 0
        self.energy = None

    def set_bfield(self,bfield):
        import bfield as bf
        if any(bfield != self.bfield):
            self.bfield = bfield
            self.H_int = [bf.calc_H_int(bfield,self.H_int_B0[n],self.xyz[n])
                          for n in range(len(self.H_int_B0))]
            self.H_hop = [bf.calc_H_hop(bfield,self.H_hop_B0[n],self.xyz[n],self.xyz[n+1])
                          for n in range(len(self.H_hop_B0))]

    def set_disorder(self,disorder,seed=None):
        import time
        self.disorder = disorder
        if seed is None:
            import time
            self.disorder_seed = long(time.time() * 256)
        else:
            self.disorder_seed = seed

    def set_energy(self,energy):
        self.energy = energy

    def H_full(self):
        res = Matrix(zeros((sum(self.N),sum(self.N)),'D'))
        N = self.N[0]
        N_old = 0
        res[0:N,0:N] = self.H_int[0]
        for i in range(1,len(self.H_int)):
            N_veryold = N_old
            N_old = N
            N += self.N[i]
            res[N_veryold:N_old,N_old:N] = self.H_hop[i-1]
            res[N_old:N,N_veryold:N_old] = adj(self.H_hop[i-1])
            res[N_old:N,N_old:N] = self.H_int[i]
        if self.disorder != 0:
            import random
            rng = random.Random(self.disorder_seed)
            for n in range(N):
                res[n,n] += (rng.random()-0.5) * self.disorder
        return res

    def eigvals(self):
        return sort(real(linalg.eigvals(self.H_full())))

    def E_H_eff(self):
        if self.disorder != 0:
            import random
            rng = random.Random(self.disorder_seed)
        assert len(self.N) > 1
        E_H_eff_LL = self.energy * Matrix(eye(self.N[0])) - self.H_int[0]
        E_H_eff_LR = - self.H_hop[0]
        E_H_eff_RR = self.energy * Matrix(eye(self.N[1])) - self.H_int[1]
        E_H_eff_RL = - adj(self.H_hop[0])
        if self.disorder != 0:
            for i in range(self.N[0]):
                E_H_eff_LL[i,i] -= (rng.random()-0.5) * self.disorder
            for i in range(self.N[1]):
                E_H_eff_RR[i,i] -= (rng.random()-0.5) * self.disorder
        for n in range(2,len(self.N)):
            E_H_eff_LX = E_H_eff_LR
            E_H_eff_XX = E_H_eff_RR
            E_H_eff_XL = E_H_eff_RL
            E_H_eff_XR = - self.H_hop[n-1]
            E_H_eff_RR = self.energy * Matrix(eye(self.N[n])) - self.H_int[n]
            E_H_eff_RX = - adj(self.H_hop[n-1])
            if self.disorder != 0:
                for i in range(self.N[n]):
                    E_H_eff_RR[i,i] -= (rng.random()-0.5) * self.disorder
            G_XX = inv(E_H_eff_XX)
            dot_XL = G_XX*E_H_eff_XL
            dot_XR = G_XX*E_H_eff_XR
            E_H_eff_LL = E_H_eff_LL - E_H_eff_LX*dot_XL
            E_H_eff_LR = E_H_eff_LX*dot_XR
            E_H_eff_RL = E_H_eff_RX*dot_XL
            E_H_eff_RR = E_H_eff_RR - E_H_eff_RX*dot_XR
        return (E_H_eff_LL,E_H_eff_LR,E_H_eff_RL,E_H_eff_RR)

    def transmission_old(self,Sigma_L,Sigma_R):
        Gamma_L = 1j*(Sigma_L - adj(Sigma_L))
        Gamma_R = 1j*(Sigma_R - adj(Sigma_R))
        (E_H_eff_LL,E_H_eff_LR,E_H_eff_RL,E_H_eff_RR) = self.E_H_eff()
        E_H_eff = r_[
            c_[E_H_eff_LL - Sigma_L,E_H_eff_LR],
            c_[E_H_eff_RL,E_H_eff_RR - Sigma_R]
        ]
        Gc = inv(E_H_eff)
        Gc_LR = Gc[0:len(Sigma_L),len(Sigma_L):len(Gc)]
        transmission = real(trace(Gamma_L*Gc_LR*Gamma_R*adj(Gc_LR)))
        return transmission

    def transmission_new(self,Sigma_0_L,Sigma_0_R,Sigma_0_L_offdiag=None,Sigma_0_R_offdiag=None):
        if self.disorder != 0:
            import random
            rng = random.Random(self.disorder_seed)
        E_H_int = []
        E_H_hop = []
        for n in range(len(self.N)):
            E_H_int.append(self.energy*Matrix(eye(self.N[n])) - self.H_int[n])
            if n>0:
                E_H_hop.append(- self.H_hop[n-1])
        disorder_range = range(len(self.N))
        if param.NO_DISORDER_IN_CONTACTS:
            disorder_range = range(len(Sigma_0_L),len(self.N)-len(Sigma_0_R))
        for n in disorder_range:
            if self.disorder != 0:
                if param.DISORDER_TYPE == 'diagonal':
                    for i in range(self.N[n]):
                        E_H_int[n][i,i] -= (rng.random()-0.5) * self.disorder
                elif param.DISORDER_TYPE == 'hopping':
                    for i,j in zip(*(E_H_int[n].nonzero())):
                        if i<j:
                            E_H_int[n][i,j] += (rng.random()-0.5) * self.disorder
                            E_H_int[n][j,i] = conj(E_H_int[n][i,j])
                    if n>0:
                        for i,j in zip(*(E_H_hop[n-1].nonzero())):
                            E_H_hop[n-1][i,j] += (rng.random()-0.5) * self.disorder
                else:
                    raise Error, 'Unknown DISORDER_TYPE: "%s"'%param.DISORDER_TYPE

        assert(len(Sigma_0_L) + len(Sigma_0_R) <= len(self.N) + 1)
        if Sigma_0_L_offdiag is not None:
            Sigma_01_L,Sigma_10_L = Sigma_0_L_offdiag
            assert len(Sigma_01_L) == len(Sigma_10_L)
            assert len(Sigma_01_L) == len(Sigma_0_L)-1
        if Sigma_0_R_offdiag is not None:
            Sigma_01_R,Sigma_10_R = Sigma_0_R_offdiag
            assert len(Sigma_01_R) == len(Sigma_10_R)
            assert len(Sigma_01_R) == len(Sigma_0_R)-1

        FINAL_CELL = (len(Sigma_0_L) + len(self.N) - len(Sigma_0_R))/2

        Sigma_L = Sigma_0_L[0]
        for n in range(len(Sigma_0_L)-1):
            E_H_eff_L = E_H_int[n] - Sigma_L
            Gs_L = inv(E_H_eff_L)
            if Sigma_0_L_offdiag is not None:
                E_H_hop_eff_01 = E_H_hop[n] - Sigma_01_L[n]
                E_H_hop_eff_10 = adj(E_H_hop[n]) - Sigma_10_L[n]
            else:
                E_H_hop_eff_01 = E_H_hop[n]
                E_H_hop_eff_10 = adj(E_H_hop[n])
            Sigma_L = E_H_hop_eff_10*Gs_L*E_H_hop_eff_01 + Sigma_0_L[n + 1]
        for n in range(len(Sigma_0_L)-1,FINAL_CELL):
            E_H_eff_L = E_H_int[n] - Sigma_L
            Gs_L = inv(E_H_eff_L)
            Sigma_L = adj(E_H_hop[n])*Gs_L*E_H_hop[n]

        Sigma_R = Sigma_0_R[-1]
        for n in range(-1,-len(Sigma_0_R),-1):
            E_H_eff_R = E_H_int[n] - Sigma_R
            Gs_R = inv(E_H_eff_R)
            if Sigma_0_R_offdiag is not None:
                E_H_hop_eff_01 = E_H_hop[n] - Sigma_01_R[n]
                E_H_hop_eff_10 = adj(E_H_hop[n]) - Sigma_10_R[n]
            else:
                E_H_hop_eff_01 = E_H_hop[n]
                E_H_hop_eff_10 = adj(E_H_hop[n])
            Sigma_R = E_H_hop_eff_01*Gs_R*E_H_hop_eff_10 + Sigma_0_R[n - 1]
        for n in range(-len(Sigma_0_R),FINAL_CELL - len(self.N),-1):
            E_H_eff_R = E_H_int[n] - Sigma_R
            Gs_R = inv(E_H_eff_R)
            Sigma_R = E_H_hop[n]*Gs_R*adj(E_H_hop[n])

        E_H_eff_C = E_H_int[FINAL_CELL] - Sigma_L - Sigma_R
        G_eff_C = inv(E_H_eff_C)
        Gamma_L = 1j*(Sigma_L - adj(Sigma_L))
        Gamma_R = 1j*(Sigma_R - adj(Sigma_R))
        T = Gamma_L*G_eff_C*Gamma_R*adj(G_eff_C)
        if param.DO_CALC_CHANNELS:
            transmission = sort(real(eigvals(T)))
            transmission = transmission[:-min(self.N)-1:-1]
        else:
            transmission = real(trace(T))
        return transmission



    def LDOS(self,Sigma_0_L,Sigma_0_R,Sigma_0_L_offdiag=None,Sigma_0_R_offdiag=None):
        if self.disorder != 0:
            import random
            rng = random.Random(self.disorder_seed)
        E_H_00 = []
        E_H_01 = []
        E_H_10 = []
        for n in range(len(self.N)):
            E_H_00.append(self.energy*Matrix(eye(self.N[n])) - self.H_int[n])
            if n>0:
                E_H_01.append(- self.H_hop[n-1])
                E_H_10.append(- adj(self.H_hop[n-1]))
        disorder_range = range(len(self.N))
        if "NO_DISORDER_IN_CONTACTS" in param:
            disorder_range = range(len(Sigma_0_L),len(self.N)-len(Sigma_0_R))
        for n in disorder_range:
            if self.disorder != 0:
                if param.DISORDER_TYPE == 'diagonal':
                    for i in range(self.N[n]):
                        E_H_00[n][i,i] -= (rng.random()-0.5) * self.disorder
                elif param.DISORDER_TYPE == 'hopping':
                    for i,j in zip(*(E_H_int[n].nonzero())):
                        if i<j:
                            E_H_00[n][i,j] += (rng.random()-0.5) * self.disorder
                            E_H_00[n][j,i] = conj(E_H_00[n][i,j])
                    if n>0:
                        for i,j in zip(*(E_H_01[n-1].nonzero())):
                            E_H_01[n-1][i,j] += (rng.random()-0.5) * self.disorder
                            E_H_10[n-1][j,i] = conj(E_H_01[n-1][i,j])
                else:
                    raise Error, 'Unknown DISORDER_TYPE: "%s"'%param.DISORDER_TYPE

        assert(len(Sigma_0_L) + len(Sigma_0_R) <= len(self.N) + 1)
        if Sigma_0_L_offdiag is not None:
            Sigma_01_L,Sigma_10_L = Sigma_0_L_offdiag
            assert len(Sigma_01_L) == len(Sigma_10_L)
            assert len(Sigma_01_L) == len(Sigma_0_L)-1
        if Sigma_0_R_offdiag is not None:
            Sigma_01_R,Sigma_10_R = Sigma_0_R_offdiag
            assert len(Sigma_01_R) == len(Sigma_10_R)
            assert len(Sigma_01_R) == len(Sigma_0_R)-1

        for n in range(len(Sigma_0_L)):
            E_H_00[n] -= Sigma_0_L[n]
            if Sigma_0_L_offdiag is not None:
                E_H_01[n] -= Sigma_01_L[n]
                E_H_10[n] -= Sigma_10_L[n]

        Sigma_R = Sigma_0_R[-1]
        for n in range(-1,-1-len(Sigma_0_R),-1):
            E_H_00[n] -= Sigma_0_R[n]
            if Sigma_0_R_offdiag is not None:
                E_H_01[n] -= Sigma_01_R[n]
                E_H_10[n] -= Sigma_10_R[n]

        res = []
        for i in range(len(self.N)):
            Sigma_L = zeros((self.N[0],)*2,'D')
            for n in range(i):
                E_H_eff_L = E_H_00[n] - Sigma_L
                Gs_L = inv(E_H_eff_L)
                Sigma_L = E_H_10[n] * Gs_L * E_H_01[n]
            Sigma_R = zeros((self.N[0],)*2,'D')
            for n in range(-1,i-len(self.N),-1):
                E_H_eff_R = E_H_00[n] - Sigma_R
                Gs_R = inv(E_H_eff_R)
                Sigma_R = E_H_01[n] * Gs_R * E_H_10[n]
            E_H_eff = E_H_00[i] - Sigma_L - Sigma_R
            G = inv(E_H_eff)
            res.append(-imag(diag(G))/pi)

        return res



def create_from_chain(chain,length,epsilon_L=0,epsilon_R=0):
    xyz = None
    if hasattr(chain,'xyz'):
	xyz = [ chain.xyz ]
    E = eye(len(chain.H_B0[0]))
    H_int = [ chain.H_B0[0]+E*epsilon_L ]
    H_hop = []
    for i in range(1,length):
	if hasattr(chain,'xyz'):
    	    xyz.append(chain.xyz.shift(i*chain.xyz.period))
        H_int.append(chain.H_B0[0]+E*(epsilon_L+i*(epsilon_R-epsilon_L)/(length-1.0)))
        H_hop.append(chain.H_B0[1])
    return conductor(xyz,H_int,H_hop)

def create_from_aperiodic_pi_orbital(aperiodic):
    xyz = [ c for c in aperiodic.cells ]
    CC_DIST = param.GRAPHENE_CC_DISTANCE
    NN_HOP = param.GRAPHENE_1STNN_HOPPING

    def hopping(pos_a,pos_b):
        if norm(pos_a - pos_b) < CC_DIST*1.35:
            return -NN_HOP
        return 0.0

#    Z_CUTOFF = CC_DIST*1.5

    H_int = []
    for i in range(len(xyz)):
        at = xyz[i].atoms
        N = len(at)
        h_int = Matrix(zeros((N,N),'D'))
        for m in range(N):
          for n in range(m+1,N):
#            if at[n].pos[2] - at[m].pos[2] > Z_CUTOFF:
#                break
            hop = hopping(at[m].pos,at[n].pos)
            h_int[m,n] = hop
            h_int[n,m] = conj(hop)
        H_int.append(h_int)
    H_hop = []
    for i in range(1,len(xyz)):
        at_A = xyz[i-1].atoms
        at_B = xyz[i].atoms
        M = len(at_A)
        N = len(at_B)
        h_hop = Matrix(zeros((M,N),'D'))
        for m in range(M):
          for n in range(N):
#            if at_B[n].pos[2] - at_A[m].pos[2] > Z_CUTOFF:
#                break
            hop = hopping(at_A[m].pos,at_B[n].pos)
            h_hop[m,n] = hop
        H_hop.append(h_hop)

    return conductor(xyz,H_int,H_hop)

def create_from_aperiodic_triozon(aperiodic):
    # based on the parametrization described in
    # doi:10.1103/PhysRevB.64.121401

    xyz = [ c for c in aperiodic.cells ]
    CC_DIST = param.GRAPHENE_CC_DISTANCE
    NN_HOP = param.GRAPHENE_1STNN_HOPPING
    TRIO_CUTOFF = param.TRIOZON_CUTOFF
    Z_CUTOFF = param.TRIOZON_Z_CUTOFF

    def hopping(pos_a,pos_b):
        if abs(norm(pos_a[:2])-norm(pos_b[:2])) < CC_DIST*0.1:
            if norm(pos_a - pos_b) < CC_DIST*1.1:
                return -NN_HOP
        else:
            d = norm(pos_b-pos_a);
            if d < TRIO_CUTOFF:
                cos_theta = vdot(pos_a[:2],pos_b[:2])/(norm(pos_a[:2])*norm(pos_b[:2]));
                return -param.TRIOZON_BETA * cos_theta * exp((param.TRIOZON_A - d)/param.TRIOZON_DELTA);
        return 0.0


    H_int = []
    for i in range(len(xyz)):
        at = xyz[i].atoms
        N = len(at)
        h_int = Matrix(zeros((N,N),'D'))
        for m in range(N):
          for n in range(m+1,N):
            if at[n].pos[2] - at[m].pos[2] > Z_CUTOFF:
                break
            hop = hopping(at[m].pos,at[n].pos)
            h_int[m,n] = hop
            h_int[n,m] = conj(hop)
        H_int.append(h_int)
    H_hop = []
    for i in range(1,len(xyz)):
        at_A = xyz[i-1].atoms
        at_B = xyz[i].atoms
        M = len(at_A)
        N = len(at_B)
        h_hop = Matrix(zeros((M,N),'D'))
        for m in range(M):
          for n in range(N):
            if at_B[n].pos[2] - at_A[m].pos[2] > Z_CUTOFF:
                break
            hop = hopping(at_A[m].pos,at_B[n].pos)
            h_hop[m,n] = hop
        H_hop.append(h_hop)

    return conductor(xyz,H_int,H_hop)


if __name__ == "__main__":
    import xyz
    import chain

    param.setdefaults()
    x = xyz.square_ladder(2)
    ch = chain.chain(x)
    cd = conductor(ch,length=3)
    cd.set_energy(1*eV)
    print cd.transmission_new(1j*Matrix(eye(2)),1j*Matrix(eye(2)))
