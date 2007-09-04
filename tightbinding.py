#!/usr/bin/env python

import re
from itertools import count

import xyz, cnt, chain, sheet
from param import param
from units import *
from calc import *

###########################

param.createdefault("GRAPHENE_1STNN_HOPPING", 2.66*eV)

def tight_binding_1stNN_graphene(xyz_coords,do_cache=True):
    maxdist = param.GRAPHENE_CC_DISTANCE * 1.1
    gamma = param.GRAPHENE_1STNN_HOPPING

    at = xyz_coords.atoms
    period = xyz_coords.period
    N = len(at)

    if isinstance(xyz_coords,xyz.chain):
        H = [ matrix(zeros((N,N))) for i in range(2) ]

        for i in range(N):
            for j in range(i+1,N):
                if norm(at[i].pos - at[j].pos) < maxdist:
                    H[0][i,j] = -gamma
                    H[0][j,i] = -gamma

        for i in range(N):
            for j in range(N):
                if norm(at[i].pos - (at[j].pos + period)) < maxdist:
                    H[1][i,j] = -gamma

        return chain.chain(H,xyz_coords,do_cache=do_cache)

    elif isinstance(xyz_coords,xyz.sheet):
        H = {}
        H[0,0] = matrix(zeros((N,N)))

        for i in range(N):
            for j in range(i+1,N):
                if norm(at[i].pos - at[j].pos) < maxdist:
                    H[0,0][i,j] = -gamma
                    H[0,0][j,i] = -gamma

        for i0,i1 in [(0,1),(1,1),(1,0),(1,-1)]:
            shift = i0 * period[0] + i1 * period[1]

            H_hop = matrix(zeros((N,N)))
            nonzero = False

            for i in range(N):
                for j in range(N):
                    if norm(at[i].pos - (at[j].pos + shift)) < maxdist:
                        H_hop[i,j] = -gamma
                        nonzero = True
            if nonzero:
                H[i0,i1] = H_hop

        return sheet.sheet(H,xyz_coords,do_cache=do_cache)



###########################

param.createdefault("TRIOZON_BETA", param.GRAPHENE_1STNN_HOPPING / 8)
param.createdefault("TRIOZON_A", 0.334 * nm)
param.createdefault("TRIOZON_DELTA", 0.045 * nm)
param.createdefault("TRIOZON_CUTOFF", param.TRIOZON_A+5*param.TRIOZON_DELTA)
param.createdefault("TRIOZON_Z_CUTOFF", (param.TRIOZON_CUTOFF**2 - (0.95*param.GRAPHITE_INTERLAYER_DISTANCE)**2)**0.5)

def tight_binding_triozon(xyz_coords,do_cache=True,graphite=False):
    # based on the parametrization described in
    # doi:10.1103/PhysRevB.64.121401

    CC_DIST = param.GRAPHENE_CC_DISTANCE
    NN_HOP = param.GRAPHENE_1STNN_HOPPING

    TRIO_CUTOFF = param.TRIOZON_CUTOFF
    Z_CUTOFF = param.TRIOZON_Z_CUTOFF
    BETA = param.TRIOZON_BETA
    A = param.TRIOZON_A
    DELTA = param.TRIOZON_DELTA


    at = xyz_coords.atoms
    period = xyz_coords.period

    Natoms = len(at)

    if isinstance(xyz_coords,xyz.chain):
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

    elif isinstance(xyz_coords,xyz.sheet):
        def hopping(pos_a,pos_b):
            if abs(pos_a[2] - pos_b[2]) < CC_DIST*0.1:
                if norm(pos_a - pos_b) < CC_DIST*1.1:
                    return -NN_HOP
            else:
                d = norm(pos_b-pos_a);
                if d < TRIO_CUTOFF:
                    return -BETA * exp((A - d)/DELTA)
            return 0.0

    if isinstance(xyz_coords,xyz.chain):
        H = [ matrix(zeros((Natoms,Natoms))) ]

        for i in range(Natoms):
            for j in range(i+1,Natoms):
                hop = hopping(at[i].pos,at[j].pos)
                if hop != 0.0:
                    H[0][i,j] = hop
                    H[0][j,i] = conj(hop)

        for n in range(1,100):
            h = Matrix(zeros((Natoms,Natoms)))
            nonzero = False
            for i in range(Natoms):
                for j in range(Natoms):
                    hop = hopping(at[i].pos,at[j].pos + period*n)
                    if hop != 0.0:
                        nonzero = True
                        h[i,j] = hop
            if not nonzero:
                break
            H.append(h)

        assert n < 99

        return chain.chain(H,xyz_coords,do_cache=do_cache)

    elif isinstance(xyz_coords,xyz.sheet):

        H = {}
        H[0,0] = matrix(zeros((Natoms,Natoms)))

        for i in range(Natoms):
            for j in range(i+1,Natoms):
                hop = hopping(at[i].pos,at[j].pos)
                if hop != 0.0:
                    H[0,0][i,j] = hop
                    H[0,0][j,i] = conj(hop)

        for i0 in range(6):
            for i1 in range(-6,6):
                if i0 == 0 and i1 <= 0:
                    continue
                shift = i0 * period[0] + i1 * period[1]

                H_hop = matrix(zeros((Natoms,Natoms)))
                nonzero = False

                for i in range(Natoms):
                    for j in range(Natoms):
                        hop = hopping(at[i].pos,at[j].pos + shift)
                        if hop != 0.0:
                            H_hop[i,j] = hop
                            nonzero = True
                if nonzero:
                    assert norm(shift) <= Z_CUTOFF + norm(period[0]) + norm(period[1])
                    H[i0,i1] = H_hop

        return sheet.sheet(H,xyz_coords,do_cache=do_cache)

###########################

class papaconstantopoulos:
    def __init__(self,fname):
        f = open(fname,'r')

        def readval(name=None):
            line = f.readline()
            line = re.sub(r'\{(.*) (.*)\}',r'{\1-\2}',line)
            vals = line.split()
            if name is not None:
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
        if R>(self.rcut+5*self.screenl):
            return 0.0
        return 1./(1+exp((R-self.rcut)/self.screenl))


    def calc_pair_SK(self,R):
        if R > (self.rcut+5*self.screenl):
            return None, None, None

        cutoff = self.cutoff(R)

        drho = exp(-self.lambda_**2*R)*cutoff

        h_sss = (self.a_sss + self.b_sss*R + self.c_sss*R**2)*exp(-self.d_sss**2*R)*cutoff
        h_sps = (self.a_sps + self.b_sps*R + self.c_sps*R**2)*exp(-self.d_sps**2*R)*cutoff
        h_pps = (self.a_pps + self.b_pps*R + self.c_pps*R**2)*exp(-self.d_pps**2*R)*cutoff
        h_ppp = (self.a_ppp + self.b_ppp*R + self.c_ppp*R**2)*exp(-self.d_ppp**2*R)*cutoff

        h_sk = matrix([
            [h_sss,-h_sps, 0.0 , 0.0 ],
            [h_sps, h_pps, 0.0 , 0.0 ],
            [ 0.0 ,  0.0 ,h_ppp, 0.0 ],
            [ 0.0 ,  0.0 , 0.0 ,h_ppp],
        ])

        s_sss = (1 + self.p_sss*R + self.q_sss*R**2 + self.r_sss*R**3)*exp(-self.s_sss**2*R)*cutoff
        s_sps = (self.p_sps*R + self.q_sps*R**2 + self.r_sps*R**3)*exp(-self.s_sps**2*R)*cutoff
        s_pps = (1 + self.p_pps*R + self.q_pps*R**2 + self.r_pps*R**3)*exp(-self.s_pps**2*R)*cutoff
        s_ppp = (1 + self.p_ppp*R + self.q_ppp*R**2 + self.r_ppp*R**3)*exp(-self.s_ppp**2*R)*cutoff

        s_sk = matrix([
            [s_sss,-s_sps, 0.0 , 0.0 ],
            [s_sps, s_pps, 0.0 , 0.0 ],
            [ 0.0 ,  0.0 ,s_ppp, 0.0 ],
            [ 0.0 ,  0.0 , 0.0 ,s_ppp],
        ])
	
        return drho, h_sk, s_sk


    def calc_pair_HS(self,Rvec):
        R = norm(Rvec)

        drho, h_sk, s_sk = self.calc_pair_SK(R)
        if drho is None:
            return None, None, None

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
        U = matrix(U)

        U4 = matrix(eye(4))
        U4[1:4,1:4] = U

        return drho, U4 * h_sk * U4.T, U4 * s_sk * U4.T


    def setup_chain(self,xyz_chain,do_cache=True):
        at = xyz_chain.atoms
        for a in at:
            a.rot4 = matrix(eye(4))
            a.rot4[1:4,1:4] = a.rot

        period = xyz_chain.period
        Natoms = len(at)

        rho = zeros((Natoms,))
        H = []
        S = []
        for n in range(20):
            at_sh = xyz_chain.shift(period*n).atoms

            h = matrix(zeros((4*Natoms,4*Natoms)))
            s = matrix(zeros((4*Natoms,4*Natoms)))
            nonzero = False
            for i in range(Natoms):
                for j in range(Natoms):
                    if n==0 and i>=j:
                        continue

                    Rvec = at[i].pos - at_sh[j].pos

                    drho, h_xyz, s_xyz = self.calc_pair_HS(Rvec)
                    if drho is None:
                        continue
                    nonzero = True

                    rho[i] += drho
                    rho[j] += drho

                    h[4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * h_xyz * at[j].rot4
                    s[4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * s_xyz * at[j].rot4

            if not nonzero:
                break
            H.append(h)
            S.append(s)

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
                H[0][4*i:4*i+4,4*j:4*j+4] = H[0][4*j:4*(j+1),4*i:4*(i+1)].H
                S[0][4*i:4*i+4,4*j:4*j+4] = S[0][4*j:4*(j+1),4*i:4*(i+1)].H

        return chain.chain(H,S=S,xyz_chain=xyz_chain,do_cache=do_cache)


    def setup_sheet(self,xyz_sheet,do_cache=True):
        at = xyz_sheet.atoms
        for a in at:
            a.rot4 = matrix(eye(4))
            a.rot4[1:4,1:4] = a.rot

        period = xyz_sheet.period
        Natoms = len(at)

        H = {}
        S = {}
        rho = zeros((Natoms,))

        H[0,0] = matrix(zeros((4*Natoms,4*Natoms)))
        S[0,0] = matrix(zeros((4*Natoms,4*Natoms)))
        for i in range(Natoms):
            for j in range(i+1,Natoms):
                Rvec = at[i].pos - at[j].pos

                drho, h_xyz, s_xyz = self.calc_pair_HS(Rvec)
                if drho is None:
                    continue

                rho[i] += drho
                rho[j] += drho

                H[0,0][4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * h_xyz * at[j].rot4
                S[0,0][4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * s_xyz * at[j].rot4

        def calc_H_S_i0_i1(i0,i1):
            shift = i0 * period[0] + i1 * period[1]
            h_hop = matrix(zeros((4*Natoms,4*Natoms)))
            s_hop = matrix(zeros((4*Natoms,4*Natoms)))
            nonzero = False

            for i in range(Natoms):
                for j in range(Natoms):
                    Rvec = at[i].pos - (at[j].pos + shift)
                    drho, h_xyz, s_xyz = self.calc_pair_HS(Rvec)

                    if drho is None:
                        continue
                    nonzero = True

                    h_hop[4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * h_xyz * at[j].rot4
                    s_hop[4*i:4*i+4,4*j:4*j+4] = at[i].rot4.T * s_xyz * at[j].rot4
            if nonzero:
                H[i0,i1] = h_hop
                S[i0,i1] = s_hop
            else:
		pass

            return nonzero

        i0 = 0
        for i1 in count(1):
            nonzero = calc_H_S_i0_i1(i0,i1)
            if not nonzero:
                break
        for i0 in count(1):
            i1 = 0
            nonzero = calc_H_S_i0_i1(i0,i1)
            if not nonzero:
                break

            for i1 in count(1):
                nonzero = calc_H_S_i0_i1(i0,i1)
                if not nonzero:
                    break

            for neg_i1 in count(1):
                nonzero = calc_H_S_i0_i1(i0,-neg_i1)
                if not nonzero:
                    break

        for i in range(Natoms):
            rho3 = rho[i]**(1/3.)
            h_s = self.alpha_s + self.beta_s * rho3**2 + self.gamma_s * rho3**4 + self.chi_s * rho3**6
            h_p = self.alpha_p + self.beta_p * rho3**2 + self.gamma_p * rho3**4 + self.chi_p * rho3**6
            H[0,0][4*i,4*i] = h_s
            H[0,0][4*i+1,4*i+1] = h_p
            H[0,0][4*i+2,4*i+2] = h_p
            H[0,0][4*i+3,4*i+3] = h_p

            S[0,0][4*i,4*i] = 1.0
            S[0,0][4*i+1,4*i+1] = 1.0
            S[0,0][4*i+2,4*i+2] = 1.0
            S[0,0][4*i+3,4*i+3] = 1.0

            for j in range(i):
                H[0,0][4*i:4*(i+1),4*j:4*(j+1)] = H[0,0][4*j:4*(j+1),4*i:4*(i+1)].H
                S[0,0][4*i:4*(i+1),4*j:4*(j+1)] = S[0,0][4*j:4*(j+1),4*i:4*(i+1)].H

        return sheet.sheet(H,S=S,xyz_sheet=xyz_sheet,do_cache=do_cache)

