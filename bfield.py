from units import *
from calc import *
from param import param

def calc_H_int(bfield,H_int_B0,xyz):
    assert(shape(H_int_B0)[0] == len(xyz.atoms))
    assert(shape(H_int_B0)[1] == len(xyz.atoms))

    if type(bfield) == float:
	if param.BFIELD_DIRECTION == 'perp':
    	    B_e_hbar = bfield * array((0.,1.,0.)) * electron / hbar
	elif param.BFIELD_DIRECTION == 'par':
    	    B_e_hbar = bfield * array((0.,0.,1.)) * electron / hbar
	elif param.BFIELD_DIRECTION == 'lateral':
    	    B_e_hbar = bfield * array((1.,0.,0.)) * electron / hbar
	else:
    	    raise "Error: unknown BFIELD_DIRECTION"
    else:
	assert len(bfield) == 3:
	B_e_hbar = array(bfield) * electron / hbar

    H_int = H_int_B0.copy()
    for i in range(len(xyz.atoms)):
        for j in range(i+1,len(xyz.atoms)):
            if H_int[i,j] != 0:
                pi = xyz.atoms[i].pos
                pj = xyz.atoms[j].pos
                mid_x = (pi[0] + pj[0])/2
                mid_y = (pi[1] + pj[1])/2
                diff_y = pj[1] - pi[1]
                diff_z = pj[2] - pi[2]
                A_e_hbar_y = B_e_hbar[2]*mid_x
                A_e_hbar_z = B_e_hbar[0]*mid_y - B_e_hbar[1]*mid_x
                factor = exp(1j*(A_e_hbar_y*diff_y + A_e_hbar_z*diff_z))
                H_int[i,j] = factor * H_int[i,j]
                H_int[j,i] = conj(H_int[i,j])
    return H_int

def calc_H_hop(bfield,H_hop_B0,xyz_0,xyz_1):
    assert(shape(H_hop_B0)[0] == len(xyz_0.atoms))
    assert(shape(H_hop_B0)[1] == len(xyz_1.atoms))

    if type(bfield) == float:
	if param.BFIELD_DIRECTION == 'perp':
    	    B_e_hbar = bfield * array((0.,1.,0.)) * electron / hbar
	elif param.BFIELD_DIRECTION == 'par':
    	    B_e_hbar = bfield * array((0.,0.,1.)) * electron / hbar
	elif param.BFIELD_DIRECTION == 'lateral':
    	    B_e_hbar = bfield * array((1.,0.,0.)) * electron / hbar
	else:
    	    raise "Error: unknown BFIELD_DIRECTION"
    else:
	assert len(bfield) == 3:
	B_e_hbar = array(bfield) * electron / hbar

    H_hop = H_hop_B0.copy()

    for i in range(len(xyz_0.atoms)):
        for j in range(len(xyz_1.atoms)):
            if H_hop[i,j] != 0:
                pi = xyz_0.atoms[i].pos
                pj = xyz_1.atoms[j].pos
                mid_x = (pi[0] + pj[0])/2
                mid_y = (pi[1] + pj[1])/2
                diff_y = pj[1] - pi[1]
                diff_z = pj[2] - pi[2]
                A_e_hbar_y = B_e_hbar[2]*mid_x
                A_e_hbar_z = B_e_hbar[0]*mid_y - B_e_hbar[1]*mid_x
                factor = exp(1j*(A_e_hbar_y*diff_y + A_e_hbar_z*diff_z))
                H_hop[i,j] = factor * H_hop[i,j]
    return H_hop
