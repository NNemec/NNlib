#!/usr/bin/env python

import xyz, cnt, chain, sheet
from param import param
from units import *
from calc import *

param.createdefault("GRAPHENE_1STNN_HOPPING", 2.66*eV)

def tight_binding_1stNN_graphene(xyz_coords,do_cache=True):
    N = len(at)

    maxdist = param.GRAPHENE_CC_DISTANCE * 1.1
    gamma = param.GRAPHENE_1STNN_HOPPING

    at = xyz_coords.atoms
    period = xyz_coords.period

    if isinstance(xyz_coords,xyz.chain):
        H = [ matrix(zeros((N,N),'D')) for i in range(2) ]

        for i in range(N):
            for j in range(i+1,N):
                if norm(at[i].pos - at[j].pos) < maxdist:
                    H[0][i,j] = -gamma
                    H[0][j,i] = -gamma

        for i in range(N):
            for j in range(N):
                if norm(at[i].pos - (at[j].pos + period)) < maxdist:
                    H[1][i,j] = -gamma

        return chain(H,xyz_coords,do_cache=do_cache)

    elif isinstance(xyz_coords,xyz.sheet):
        H = {}
        H[0,0] = matrix(zeros((N,N),'D'))

        for i in range(N):
            for j in range(i+1,N):
                if norm(at[i].pos - at[j].pos) < maxdist:
                    H[0,0][i,j] = -gamma
                    H[0,0][j,i] = -gamma

        for i0,i1 in [(0,1),(1,1),(1,0),(1,-1)]:
            shift = i0 * period[0] + i1 * period[1]

            H_hop = matrix(zeros((N,N),'D'))
            nonzero = False

            for i in range(N):
                for j in range(N):
                    if norm(at[i].pos - (at[j].pos + shift)) < maxdist:
                        H_hop[i,j] = -gamma
                        nonzero = True
            if nonzero:
                H[i0,i1] = H_hop

        return sheet(H,xyz_coords,do_cache=do_cache)



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
        H = [ matrix(zeros((Natoms,Natoms),'D')) ]

        for i in range(Natoms):
            for j in range(i+1,Natoms):
                hop = hopping(at[i].pos,at[j].pos)
                if hop != 0.0:
                    H[0][i,j] = hop
                    H[0][j,i] = conj(hop)

        for n in range(1,100):
            h = Matrix(zeros((Natoms,Natoms),'D'))
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
        H[0,0] = matrix(zeros((Natoms,Natoms),'D'))

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

                H_hop = matrix(zeros((Natoms,Natoms),'D'))
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

