#!/usr/bin/env python

from calc import *
from units import *
import re
from chain import chain

class papa:
    def __init__(self,fname):
        f = open(fname,'r')

        def readval(name=None):
            line = f.readline()
            line = re.sub(r'\{(.*) (.*)\}',r'{\1-\2}',line)
#           print line
            vals = line.split()
            if name is not None:
#                print vals[3],name
                assert vals[3] == name.split()[0]
            val = float(vals[0])
            return val

        for i in range(3):
            f.readline()

        rcut,screenl = f.readline().split()[:2]
        self.rcut = float(rcut) * bohr
        self.screenl = float(screenl) * bohr

        for i in range(3):
            f.readline()

        setattr(self,'lambda_',readval('lambda')/bohr**.5)

        for l in ['s','p']:
            setattr(self,'alpha_%s'%l, readval('a_%s'%l)*rydberg)
            setattr(self,'beta_%s'%l,  readval('b_%s'%l)*rydberg)
            setattr(self,'gamma_%s'%l, readval('c_%s'%l)*rydberg)
            setattr(self,'chi_%s'%l, readval('d_%s'%l)*rydberg)

        for i in range(8):
            f.readline()

        for l in ['sss','sps','pps','ppp']:
            ltrans = '{%s-%s}'%(l[:2],{'s':'sigma','p':'pi'}[l[2]])
            setattr(self,'a_%s'%l, readval('e_%s'%ltrans)*rydberg)
            setattr(self,'b_%s'%l, readval('f_%s'%ltrans)*rydberg/bohr)
            setattr(self,'c_%s'%l, readval('fbar_%s'%ltrans)*rydberg/bohr**2)
            setattr(self,'d_%s'%l, readval('g_%s'%ltrans)/bohr**.5)

        for i in range(24):
            f.readline()

        for l in ['sss','sps','pps','ppp']:
            ltrans = '{%s-%s}'%(l[:2],{'s':'sigma','p':'pi'}[l[2]])
            setattr(self,'p_%s'%l, readval('e_%s'%ltrans)/bohr)
            setattr(self,'q_%s'%l, readval('f_%s'%ltrans)/bohr**2)
            setattr(self,'r_%s'%l, readval('fbar_%s'%ltrans)/bohr**3)
            setattr(self,'s_%s'%l, readval('g_%s'%ltrans)/bohr**.5)

        for i in range(24):
            f.readline()

        assert f.readline() == ""

        f.close()

    def cutoff(self,R):
        return 1./(1+exp((R-self.rcut)/self.screenl))


    def setup_chain(self,xyz,do_cache=True):
        at = xyz.atoms
        for a in at:
            a.rot4 = Matrix(eye(4))
            a.rot4[1:4,1:4] = a.rot

        period = xyz.period

        Natoms = len(xyz.atoms)

        rho = zeros((Natoms,),'d')
        H = []
        S = []
        for n in range(20):
            at_sh = xyz.shift(period*n).atoms

            h = Matrix(zeros((4*Natoms,4*Natoms),'D'))
            s = Matrix(zeros((4*Natoms,4*Natoms),'D'))
            nonzero = False
            for i in range(Natoms):
                for j in range(Natoms):
                    if n==0 and i>=j:
                        continue

                    Rvec = at[i].pos - at_sh[j].pos

                    R = norm(Rvec)
                    if R > (self.rcut+5*self.screenl):
                        continue

                    U = zeros((3,3))

                    U[:,0] = Rvec/R
                    l=0
                    if abs(U[1,0])<abs(U[0,0]):
                        l=1
                    if abs(U[2,0])<abs(U[l,0]):
                        l=2
                    U[l,1] = 1
                    U[:,2] = cross(U[:,0],U[:,1])
                    U[:,2] = U[:,2]/norm(U[:,2])
                    U[:,1] = cross(U[:,0],U[:,2])
                    U = Matrix(U)

                    nonzero = True
                    cutoff = self.cutoff(R)
                    drho = exp(-self.lambda_**2*R)*cutoff
                    rho[i] += drho
                    rho[j] += drho

                    U4 = Matrix(eye(4))
                    U4[1:4,1:4] = U
		    

                    h_sss = (self.a_sss + self.b_sss*R + self.c_sss*R**2)*exp(-self.d_sss**2*R)*cutoff
                    h_sps = (self.a_sps + self.b_sps*R + self.c_sps*R**2)*exp(-self.d_sps**2*R)*cutoff
                    h_pps = (self.a_pps + self.b_pps*R + self.c_pps*R**2)*exp(-self.d_pps**2*R)*cutoff
                    h_ppp = (self.a_ppp + self.b_ppp*R + self.c_ppp*R**2)*exp(-self.d_ppp**2*R)*cutoff

                    h_sk = Matrix([
                        [h_sss,h_sps, 0.0 , 0.0 ],
                        [h_sps,h_pps, 0.0 , 0.0 ],
                        [ 0.0 , 0.0 ,h_ppp, 0.0 ],
                        [ 0.0 , 0.0 , 0.0 ,h_ppp],
                    ])

                    h[4*i:4*(i+1),4*j:4*(j+1)] = at[i].rot4.T * U4 * h_sk * U4.T * at[j].rot4


                    s_sss = (1 + self.p_sss*R + self.q_sss*R**2 + self.r_sss*R**3)*exp(-self.s_sss**2*R)*cutoff
                    s_sps = (self.p_sps*R + self.q_sps*R**2 + self.r_sps*R**3)*exp(-self.s_sps**2*R)*cutoff
                    s_pps = (1 + self.p_pps*R + self.q_pps*R**2 + self.r_pps*R**3)*exp(-self.s_pps**2*R)*cutoff
                    s_ppp = (1 + self.p_ppp*R + self.q_ppp*R**2 + self.r_ppp*R**3)*exp(-self.s_ppp**2*R)*cutoff

                    s_sk = Matrix([
                        [s_sss,s_sps, 0.0 , 0.0 ],
                        [s_sps,s_pps, 0.0 , 0.0 ],
                        [ 0.0 , 0.0 ,s_ppp, 0.0 ],
                        [ 0.0 , 0.0 , 0.0 ,s_ppp],
                    ])

                    s[4*i:4*(i+1),4*j:4*(j+1)] = at[i].rot4.T * U4 * s_sk * U4.T * at[j].rot4

            if nonzero:
                H.append(h)
                S.append(s)
            else:
                break

        for i in range(Natoms):
            rho3 = rho[i]**(1/3.)
            h_s = self.alpha_s + self.beta_s * rho3**2 + self.gamma_s * rho3**4 + self.chi_s * rho3**6
            h_p = self.alpha_p + self.beta_p * rho3**2 + self.gamma_p * rho3**4 + self.chi_p * rho3**6
            H[0][4*i,4*i] = h_s
            H[0][4*i+1,4*i+1] = h_p
            H[0][4*i+2,4*i+2] = h_p
            H[0][4*i+3,4*i+3] = h_p

            S[0][4*i,4*i] = 1.0
            S[0][4*i+1,4*i+1] = 1.0
            S[0][4*i+2,4*i+2] = 1.0
            S[0][4*i+3,4*i+3] = 1.0

            for j in range(i):
                H[0][4*i:4*(i+1),4*j:4*(j+1)] = transpose(H[0][4*j:4*(j+1),4*i:4*(i+1)])
                S[0][4*i:4*(i+1),4*j:4*(j+1)] = transpose(S[0][4*j:4*(j+1),4*i:4*(i+1)])

        return chain(H,S=S,xyz_chain=xyz,do_cache=do_cache)
