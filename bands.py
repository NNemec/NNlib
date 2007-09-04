#!/usr/bin/env python

from scan import scan_adaptive
from numpy import *
from units import *

def floor(x):
    return int(x//1.0)

def ceil(x):
    return -int(-x//1.0)

class scan_bands(scan_adaptive):
    def __init__(self,
        chain,
        verbose = False,
	precision = 1e-3,
    ):
        scan_adaptive.__init__(self,
            chain.band_energies,
            initial_xgrid = linspace(0,2*pi,9)[:-1],
            period = 2*pi,
            verbose = verbose,
            sameyscale = True,
	    precision = precision,
        )
        self.Nbands = self.y.shape[1]

    def calc_N_bands(self,energy):
        kcut, = self.find_valuecut(
            energy,
        )
        N_bands = len(kcut)
        return N_bands


    def calc_DOS(self,energy):
        """
        returns density of states per unit cell
        """
        slopecut, = self.find_valuecut(
            energy,
            calccutx = False,
            calcslope = True,
        )
        DOS = sum(abs(1/slopecut) / (2*pi))
        return DOS


    def calc_NOS(self,maxenergy):
        """
        returns number of states per unit cell up to energy maxenergy
        """
        kcut,signcut = self.find_valuecut(
            value=maxenergy,
            calccutx=True,
            calcslopesign=True,
        )
        NOS = sum(signcut*kcut) / (2*pi)
        NOS += sum(self.y[0,:] < maxenergy) # to add bands that wrap to next brillouin zone
#        print maxenergy / eV, NOS
        return NOS

    def find_fermi_energy(self,NOE,yminstep=None):
        """
        calculate the Fermi energy for a given number of electrons per unit cell (NOS)
        """
	
	
        print "NOE:",NOE
	Emax = self.y.max()
	Emin = self.y.min()

        if yminstep is None:
            yminstep = self.yminstep
        if yminstep is None:
            yminstep = (Emax - Emin) * self.precision

        NOSmax = self.calc_NOS(Emax)
        NOSmin = self.calc_NOS(Emin)

        while (Emax - Emin) > yminstep or (NOSmax-NOSmin) > self.precision:
#            Emid = (Emax+Emin) / 2
            Emid = (NOE-NOSmin) * (Emax-Emin) / (NOSmax-NOSmin) + Emin
            NOSmid = self.calc_NOS(Emid)
	    assert Emin <= Emid 
	    assert Emid <= Emax 
            print Emax/eV,Emin/eV,NOSmax,NOSmin,Emid/eV,NOSmid
            if NOSmid > NOE:
		if Emax - Emid < yminstep:
		    return Emid
                Emax = Emid
                NOSmax = NOSmid
            elif NOSmid < NOE:
		if Emid - Emin < yminstep:
		    return Emid
                Emin = Emid
                NOSmin = NOSmid
            else:
                return Emid

        return (Emin+Emax) / 2

if __name__ == "__main__":
    import cnt, tightbinding, units
    import chain as NNlib_chain

    papa = tightbinding.papaconstantopoulos("c_par.105.tbparam")

    xyz = cnt.GNR_zigzag(6)
    chain = papa.setup_chain(xyz)
    chain = NNlib_chain.chain([h[2::4,2::4] for h in chain.H],chain.xyz,S=[s[2::4,2::4] for s in chain.S])
    chain = chain.multiply()
    
    scan = scan_bands(
	chain,
	precision=1e-3,
	verbose=True,
    )
    scan.do_scan()
    scan.sort_crossing()

    band_k = scan.x
    band_energy = scan.y


    if scan.period is not None:
        band_k = concatenate((
            band_k,
            band_k[:1] + scan.period,
        ))
	    
        addyvals = band_energy[:1,:]
        if hasattr(scan,'totalpermut'):
            addyvals = addyvals[:,scan.totalpermut.argsort()]
        band_energy = concatenate((
            band_energy,
            addyvals,
        ),axis=0)

    Emin = band_energy.min()
    Emax = band_energy.max()
    Emid = (Emin+Emax)*.5
    Emin = (Emin-Emid)*1.1 + Emid
    Emax = (Emax-Emid)*1.1 + Emid
    
    def calc_DOS(e):
	return scan.calc_DOS(e)
    
    scan_DOS = scan_adaptive(
	calc_DOS,
	linspace(Emin,Emax,5),
	verbose = True,
	precision = 1e-3,
    )
    scan_DOS.do_scan()

    def calc_NOS(e):
	return scan.calc_NOS(e)
    
    scan_NOS = scan_adaptive(
	calc_NOS,
	linspace(Emin,Emax,5),
	verbose = True,
	precision = 1e-3,
    )
    scan_NOS.do_scan()

    NOE = arange(scan.Nbands+1)
    E_F = array([scan.find_fermi_energy(n) for n in NOE])

    import pylab
    pylab.subplot(1,3,1)
    for b in range(scan.Nbands):
        pylab.plot(band_k,band_energy[:,b] / eV,'.-')
    for i in range(len(NOE)):
        print "%i: %g eV, %g electrons"%(NOE[i],E_F[i]/units.eV,scan.calc_NOS(E_F[i]))
        pylab.axhline(E_F[i] / eV)

    pylab.subplot(1,3,2)
    pylab.plot(scan_DOS.y,scan_DOS.x / eV,'.')
    for i in range(10):
	scan_DOS.reduce_visible()
    pylab.plot(scan_DOS.y,scan_DOS.x / eV,'.-')

    pylab.subplot(1,3,3)
    pylab.plot(scan_NOS.y,scan_NOS.x / eV,'.')
    for i in range(10):
        scan_NOS.reduce_visible()
    pylab.plot(scan_NOS.y,scan_NOS.x / eV,'.-')

    pylab.show()
