from units import *
from calc import *

CUTOFF_REL = 7

def conductance(transmission,energy,temperature,at_energy):
    if temperature == 0:
	if at_energy > energy[-1]:
	    raise "error: at_energy too high"
	e = 0
	while energy[e] < at_energy:
	    e += 1
	
	if e == 0:
	    raise "error: at_energy too low"
	return (((energy[e] - at_energy)*transmission[e] + (at_energy - energy[e-1])*transmission[e-1])/(energy[e] - energy[e-1])),""

    beta = 1/(k_B*temperature)

    ENERGY_CUTOFF = CUTOFF_REL/beta
    def f(d_E):
        return ((ENERGY_CUTOFF-abs(d_E))>0)*beta/((exp(beta*d_E)+1)*(exp(-beta*d_E)+1))

    datapoints = 0
    sum = 0.
    control = 0.
    for e in range(1,len(energy)-1):
        d_E = abs(energy[e]-at_energy)
        if d_E < ENERGY_CUTOFF:
            datapoints += 1
            factor = beta/((exp(beta*d_E)+1)*(exp(-beta*d_E)+1)) * (energy[e+1]-energy[e-1])/2
            sum += transmission[e] * factor
            control += factor

    msg = "datapoints: %i, control: %g"%(datapoints,control)
    if at_energy - energy[0] < ENERGY_CUTOFF:
        msg += "overflow at minimum"

    if energy[-1] - at_energy < ENERGY_CUTOFF:
        msg += "overflow at maximum"

    conductance = sum/control*G_0

    return conductance,msg
