import tables as pytables
from units import *

class DataDummy:
    def __init__(self):
        self.__dict__['_v_attrnamesuser'] = {}

    def __getattr__(self,name):
        return self._v_attrnamesuser[name]

    def __setattr__(self,name,value):
        self._v_attrnamesuser[name] = value

    def __delattr__(self,name):
        del self._v_attrnamesuser[name]

class Param:
    def setdefaults(self):
        data = self.data

        for name in data._v_attrnamesuser:
            delattr(data,name)

        data.LOPEZ_SANCHO_ETA = 1e-5*eV
        data.LOPEZ_SANCHO_EPSILON = 1e-4*eV

        data.GRAPHENE_CC_DISTANCE = 1.4226*angstrom
        data.GRAPHITE_INTERLAYER_DISTANCE = 3.44*angstrom
        data.GRAPHENE_1STNN_HOPPING = 2.66*eV

        data.LATTICE_CONSTANT = data.GRAPHENE_CC_DISTANCE
        data.HOPPING_CONSTANT = data.GRAPHENE_1STNN_HOPPING

#        data.LEAD_TYPE = "wideband"
#        data.LEAD_TYPE = "coating_wideband"
#        data.WIDEBAND_ENERGY = 1.0*eV
        data.LEAD_TYPE = "lopez_sancho"

#        data.DO_CALC_CHANNELS = True
#        data.NO_DISORDER_IN_CONTACTS = True
        data.BFIELD_IN_LEADS = True
        data.BFIELD_DIRECTION = 'perp'
#        data.BFIELD_DIRECTION = 'par'

        data.DISORDER_TYPE = "diagonal"
#        data.DISORDER_TYPE = "hopping"

        data.TRIOZON_BETA = data.GRAPHENE_1STNN_HOPPING / 8
        data.TRIOZON_A = 0.334 * nm
        data.TRIOZON_DELTA = 0.045 * nm
        data.TRIOZON_CUTOFF = data.TRIOZON_A+5*data.TRIOZON_DELTA
        data.TRIOZON_Z_CUTOFF = (data.TRIOZON_CUTOFF**2 - (0.95*data.GRAPHITE_INTERLAYER_DISTANCE)**2)**0.5

    def reset(self):
        self.__dict__['temp'] = {}
        self.__dict__['data'] = DataDummy()

    def __init__(self):
        self.reset()

    def activate(self,data,copy=False):
        if isinstance(data,pytables.Group):
            data = data._v_attrs
        elif isinstance(data,pytables.Leaf):
            data = data.attrs
        if copy:
            for name in data._v_attrnamesuser:
                delattr(data,name)
            for name in self.data._v_attrnamesuser:
                setattr(data,name,getattr(self.data,name))
        self.__dict__['data'] = data

    def settemp(self,name,value=None):
        self.temp[name] = value

    def __getattr__(self,name):
        if name in self.temp:
            return self.temp[name]
        else:
            return getattr(self.data,name)

    def __setattr__(self,name,value):
        if name in self.temp:
            self.temp[name] = value
        else:
            setattr(self.data,name,value)

    def __delattr__(self,name):
        if name in self.temp:
            del self.temp[name]
        else:
            delattr(self.data,name)

    def __contains__(self,name):
        return name in self.temp or hasattr(self.data,name)

param = Param()
