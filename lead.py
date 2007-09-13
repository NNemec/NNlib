from calc import *
from units import *

class wideband:
    def __init__(self,contact,factor=1.0/eV):
        self.contact = contact
        s = contact.shape
        assert(s[0] == s[1])
        N = s[0]
        self._Gs = matrix(eye(N))*1.0j*factor
        self._Sigma_L = self.contact.H * self._Gs * self.contact
        self._Sigma_R = self.contact * self._Gs * self.contact.H

    def Gs_L(self,energy=None,bfield=None):
        return self._Gs

    def Gs_R(self,energy=None,bfield=None):
        return self._Gs

    def Sigma_L(self,energy=None,bfield=None):
        return self._Sigma_L

    def Sigma_R(self,energy=None,bfield=None):
        return self._Sigma_R

    def set_energy(self, energy):
        pass

    def set_bfield(self, bfield):
        pass

class nondiag_wideband_L:
    def __init__(self,contact00,contact01,contact10,SDOS=1.0/eV):
        assert len(contact00) == len(contact01)+1
        assert len(contact00) == len(contact10)+1

        self._Gs = []
        for n in range(len(contact00)):
            N = contact00[n].shape[0]
            self._Gs.append(matrix(eye(N))*1.0j*SDOS)
        self._Sigma_L_00 = []
        for n in range(len(contact00)):
            self._Sigma_L_00.append(contact00[n].H*self._Gs[n]*contact00[n])
        for n in range(1,len(contact00)):
            self._Sigma_L_00[n] += contact01[n-1].H*self._Gs[n-1]*contact01[n-1]
        for n in range(len(contact00)-1):
            self._Sigma_L_00[n] += contact10[n].H*self._Gs[n+1]*contact10[n]

        self._Sigma_L_01 = []
        for n in range(len(contact00)-1):
            self._Sigma_L_01.append(
                contact00[n].H*self._Gs[n]*contact01[n]
            +   contact10[n].H*self._Gs[n+1]*contact00[n+1]
            )

        self._Sigma_L_10 = []
        for n in range(len(contact00)-1):
            self._Sigma_L_10.append(
                contact01[n].H*self._Gs[n]*contact00[n]
            +   contact00[n+1].H*self._Gs[n+1]*contact10[n]
            )

    def Gs_L(self,energy=None,bfield=None):
        return self._Gs

    def Sigma_L_diag(self,energy=None,bfield=None):
        return self._Sigma_L_00

    def Sigma_L_offdiag(self,energy=None,bfield=None):
        return (self._Sigma_L_01,self._Sigma_L_10)

    def set_energy(self, energy):
        pass

    def set_bfield(self, bfield):
        pass


class nondiag_wideband_R:
    def __init__(self,contact00,contact01,contact10,SDOS=1.0/eV):
        assert len(contact00) == len(contact01)+1
        assert len(contact00) == len(contact10)+1

        self._Gs = []
        for n in range(len(contact00)):
            N = contact00[n].shape[1]
            self._Gs.append(matrix(eye(N))*1.0j*SDOS)
        self._Sigma_R_00 = []
        for n in range(len(contact00)):
            self._Sigma_R_00.append(contact00[n]*self._Gs[n]*contact00[n].H)
        for n in range(1,len(contact00)):
            self._Sigma_R_00[n] += contact10[n-1]*self._Gs[n-1]*contact10[n-1].H
        for n in range(len(contact00)-1):
            self._Sigma_R_00[n] += contact01[n]*self._Gs[n+1]*contact01[n].H

        self._Sigma_R_01 = []
        for n in range(len(contact00)-1):
            self._Sigma_R_01.append(
                contact00[n]*self._Gs[n]*contact10[n].H
            +   contact01[n]*self._Gs[n+1]*contact00[n+1].H
            )

        self._Sigma_R_10 = []
        for n in range(len(contact00)-1):
            self._Sigma_R_10.append(
                contact10[n]*self._Gs[n]*contact00[n].H
            +   contact00[n+1]*self._Gs[n+1]*contact01[n].H
            )

    def Gs_R(self,energy=None,bfield=None):
        return self._Gs

    def Sigma_R_diag(self,energy=None,bfield=None):
        return self._Sigma_R_00

    def Sigma_R_offdiag(self,energy=None,bfield=None):
        return (self._Sigma_R_01,self._Sigma_R_10)

    def set_energy(self, energy):
        pass

    def set_bfield(self, bfield):
        pass


class lopez_sancho:
    def __init__(self,chain,tunneling = 1.0,stoner_shift = 0.0*eV):
        self.chain = chain
        self.tunneling = tunneling
        self.stoner_shift = stoner_shift
        self.energy = None

    def Gs_L(self,energy=None):
        if energy is not None:
            set_energy(energy)
        return self.chain.Gs_L(self.energy-self.stoner_shift)

    def Gs_R(self,energy=None):
        if energy is not None:
            set_energy(energy)
        return self.chain.Gs_R(self.energy-self.stoner_shift)

    def Sigma_L(self,energy=None):
        return (
            self.chain.H[1].H
            * self.Gs_L(energy)
            * self.chain.H[1]
            * self.tunneling**2
        )

    def Sigma_R(self,energy=None):
        return (
            self.chain.H[1]
            * self.Gs_R(energy)
            * self.chain.H[1].H
            * self.tunneling**2
        )

    def set_energy(self, energy):
        self.energy = energy

    def set_bfield(self, bfield):
        self.chain.set_bfield(bfield)
