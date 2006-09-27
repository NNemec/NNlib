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
#	    print line
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

        setattr(self,'lambda_',readval('lambda'))

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
            setattr(self,'p_%s'%l, readval('e_%s'%ltrans)*rydberg)
            setattr(self,'q_%s'%l, readval('f_%s'%ltrans)*rydberg/bohr)
            setattr(self,'r_%s'%l, readval('fbar_%s'%ltrans)*rydberg/bohr**2)
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
        for n in range(20):
            at_sh = xyz.shift(period*n).atoms

            h_hop = Matrix(zeros((4*Natoms,4*Natoms),'D'))
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
		    i=0
		    if abs(U[1,0])<abs(U[0,0]):
			i=1
		    if abs(U[2,0])<abs(U[i,0]):
			i=2
		    U[i,1]
		    U[:,2] = cross(U[:,0],U[:,1])
		    U[:,2] = U[:,2]/norm(U[:,2])
		    U[:,1] = cross(U[:,0],U[:,2])
		    U = Matrix(U)

                    nonzero = True
                    cutoff = self.cutoff(R)
                    drho = exp(-self.lambda_**2*R)*cutoff
                    rho[i] += drho
                    rho[j] += drho

                    h_sss = (self.a_sss + self.b_sss*R + self.c_sss*R**2)*exp(-self.d_sss**2*R)*cutoff
                    h_sps = (self.a_sps + self.b_sps*R + self.c_sps*R**2)*exp(-self.d_sps**2*R)*cutoff
                    h_pps = (self.a_pps + self.b_pps*R + self.c_pps*R**2)*exp(-self.d_pps**2*R)*cutoff
                    h_ppp = (self.a_ppp + self.b_ppp*R + self.c_ppp*R**2)*exp(-self.d_ppp**2*R)*cutoff

                    sk_hop = Matrix([
                        [h_sss,h_sps, 0.0 , 0.0 ],
                        [h_sps,h_pps, 0.0 , 0.0 ],
                        [ 0.0 , 0.0 ,h_ppp, 0.0 ],
                        [ 0.0 , 0.0 , 0.0 ,h_ppp],
                    ])
		    
		    U4 = Matrix(eye(4))
		    U4[1:4,1:4] = U
		    
                    ##### sk_hop is not yet correctly oriented !! #####

                    h_hop[4*i:4*(i+1),4*j:4*(j+1)] = at[i].rot4 * U4 * sk_hop * U4.T * at[j].rot4.T
            if nonzero:
                H.append(h_hop)
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

            for j in range(i):
                H[0][4*i:4*(i+1),4*j:4*(j+1)] = transpose(H[0][4*j:4*(j+1),4*i:4*(i+1)])

        return chain(H,xyz,do_cache=do_cache)
