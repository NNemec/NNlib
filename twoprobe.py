from calc import *
from units import *
from param import param

import lead as _lead
import conductor as _conductor

class twoprobe:
    def __init__(self,conductor,lead_L,lead_R):
        self.energy = None
        self.bfield = None
        self.conductor = conductor
        assert type(lead_L) is list
        self.lead_L = lead_L
        assert type(lead_R) is list
        self.lead_R = lead_R
        self.Sigma_L = None
        self.Sigma_R = None

    def set_energy(self,energy):
        if self.energy != energy:
            self.energy = energy
            for l in self.lead_L:
                l.set_energy(energy)
            for l in self.lead_R:
                l.set_energy(energy)
            self.conductor.set_energy(energy)
            self.Sigma_L = None
            self.Sigma_R = None

    def set_bfield(self,bfield):
        if self.bfield is None or bfield != self.bfield:
            self.bfield = bfield
            self.conductor.set_bfield(bfield)
            if 'BFIELD_IN_LEADS' in param:
                for l in self.lead_L:
                    l.set_bfield(bfield)
                for l in self.lead_R:
                    l.set_bfield(bfield)
                self.Sigma_L = None
                self.Sigma_R = None

    def set_disorder(self,disorder,seed=None):
        self.conductor.set_disorder(disorder,seed)

    def transmission(self):
        if self.Sigma_L is None:
            self.Sigma_L = []
            for l in self.lead_L:
                self.Sigma_L.append(l.Sigma_L())
            self.Sigma_R = []
            for l in self.lead_R:
                self.Sigma_R.append(l.Sigma_R())
        return self.conductor.transmission_new(self.Sigma_L,self.Sigma_R)

    def LDOS(self,energy=None):
        if self.Sigma_L is None:
            self.Sigma_L = []
            for l in self.lead_L:
                self.Sigma_L.append(l.Sigma_L())
            self.Sigma_R = []
            for l in self.lead_R:
                self.Sigma_R.append(l.Sigma_R())
        return self.conductor.LDOS(self.Sigma_L,self.Sigma_R)


class twoprobe_nondiagleads:
    def __init__(self,conductor,lead_L,lead_R):
        self.energy = None
        self.bfield = None
        self.conductor = conductor
        assert isinstance(lead_L,_lead.nondiag_wideband_L)
        self.lead_L = lead_L
        assert isinstance(lead_R,_lead.nondiag_wideband_R)
        self.lead_R = lead_R
        self.Sigma_L = None
        self.Sigma_R = None

    def set_energy(self,energy):
        if self.energy != energy:
            self.energy = energy
            self.lead_L.set_energy(energy)
            self.lead_R.set_energy(energy)
            self.conductor.set_energy(energy)
            self.Sigma_L = None
            self.Sigma_R = None

    def set_bfield(self,bfield):
        if self.bfield is None or bfield != self.bfield:
            self.bfield = bfield
            self.conductor.set_bfield(bfield)
            if 'BFIELD_IN_LEADS' in param:
                self.lead_L.set_bfield(bfield)
                self.lead_R.set_bfield(bfield)
                self.Sigma_L = None
                self.Sigma_R = None

    def set_disorder(self,disorder,seed=None):
        self.conductor.set_disorder(disorder,seed)

    def transmission(self,energy=None):
        if energy is not None:
            self.set_energy(energy)
        return self.conductor.transmission_new(
            self.lead_L.Sigma_L_diag(),
            self.lead_R.Sigma_R_diag(),
            self.lead_L.Sigma_L_offdiag(),
            self.lead_R.Sigma_R_offdiag(),
        )



def create_from_chain(chain,conductor_size,contact_size_L=1,contact_size_R=1):
    conductor = _conductor.create_from_chain(chain,conductor_size)

    if param.LEAD_TYPE == 'wideband':
        lead = _lead.wideband(
            contact = where(chain.H_hop,1.0,0.0)*eV,
            factor=param.WIDEBAND_ENERGY/eV**2,
        )

    elif param.LEAD_TYPE == 'coating_wideband':
        if 'COATING_WIDEBAND_N_CONTACT' in param:
            contact = Matrix(zeros((chain.N_atoms,)*2,'D'))
            for n in range(param.COATING_WIDEBAND_N_CONTACT):
                contact[n,n] = 1.0*eV
        else:
            contact = Matrix(eye(chain.N_atoms))*1.0*eV
        lead = _lead.wideband(
            contact = contact,
            factor=param.WIDEBAND_ENERGY/eV**2,
        )
    elif param.LEAD_TYPE == 'lopez_sancho':
        tunneling = 1
        if "TUNNELING_FACTOR" in param:
            tunneling *= param.TUNNELING_FACTOR
        lead = _lead.lopez_sancho(
            chain,
            tunneling = tunneling,
        )

    return twoprobe(
        conductor,
        [lead]*contact_size_L,
        [lead]*contact_size_R,
    )


def create_from_aperiodic(aperiodic,contact_length_L,contact_length_R):
    conductor = _conductor.create_from_aperiodic_triozon(aperiodic)
    diameter = max([norm(a.pos[:2]) for a in aperiodic.cells[0].atoms])
    mindiam = diameter - param.GRAPHENE_CC_DISTANCE*0.5
    lead_L = []
    lead_R = []
    finished = False
    for n in range(len(aperiodic.cells)):
        atoms = aperiodic.cells[n].atoms
        contact = Matrix(zeros((len(atoms),)*2,'D'))
        j = 0
        for i in range(len(atoms)):
            if atoms[i].pos[2] > contact_length_L:
                finished = True
            else:
                if norm(atoms[i].pos[:2]) >= mindiam:
                    contact[j,i] = 1.0*eV
                    j += 1
        lead_L.append(_lead.wideband(
            contact[:j,:],
            factor=param.WIDEBAND_ENERGY/eV**2,
        ))
        if finished:
            break

    for n in reversed(range(len(aperiodic.cells))):
        atoms = aperiodic.cells[n].atoms
        contact = Matrix(zeros((len(atoms),)*2,'D'))
        j = 0
        for i in range(len(atoms)):
            if atoms[i].pos[2] < aperiodic.length - contact_length_R:
                finished = True
            else:
                if norm(atoms[i].pos[:2]) >= mindiam:
                    contact[i,j] = 1.0*eV
                    j += 1
        lead_R.append(_lead.wideband(
            contact[:,:j],
            factor=param.WIDEBAND_ENERGY/eV**2
        ))
        if finished:
            break
    lead_R.reverse()

    return twoprobe(
        conductor,
        lead_L,
        lead_R,
    )
