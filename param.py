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
    def reset(self):
        self.__dict__['temp'] = {}
        self.__dict__['data'] = DataDummy()

    def __init__(self):
        self.__dict__['defaults'] = {}
        self.reset()

    def setdefaults(self):
        for name in self.data._v_attrnamesuser:
            delattr(self.data,name)

        for name,value in self.defaults.iteritems():
            setattr(self.data,name,value)

    def createdefault(self,name,value):
        self.defaults[name] = value
        if not hasattr(self.data,name):
            setattr(self.data,name,value)

    def write(self,data):
        if isinstance(data,pytables.Group):
            data = data._v_attrs
        elif isinstance(data,pytables.Leaf):
            data = data.attrs
        for name in data._v_attrnamesuser:
            delattr(data,name)
        for name in self.data._v_attrnamesuser:
            setattr(data,name,getattr(self.data,name))
        self.__dict__['data'] = data

    def read(self,data):
        if isinstance(data,pytables.Group):
            data = data._v_attrs
        elif isinstance(data,pytables.Leaf):
            data = data.attrs
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


param.createdefault("GRAPHENE_CC_DISTANCE", 1.4226*angstrom)
param.createdefault("GRAPHITE_INTERLAYER_DISTANCE", 3.44*angstrom)
param.createdefault("GRAPHENE_1STNN_HOPPING", 2.66*eV)

param.createdefault("LATTICE_CONSTANT", param.GRAPHENE_CC_DISTANCE)
param.createdefault("HOPPING_CONSTANT", param.GRAPHENE_1STNN_HOPPING)

#        data.DO_CALC_CHANNELS", True
#        data.NO_DISORDER_IN_CONTACTS", True

param.createdefault("DISORDER_TYPE", "diagonal")
#        data.DISORDER_TYPE", "hopping"

param.createdefault("TRIOZON_BETA", param.GRAPHENE_1STNN_HOPPING / 8)
param.createdefault("TRIOZON_A", 0.334 * nm)
param.createdefault("TRIOZON_DELTA", 0.045 * nm)
param.createdefault("TRIOZON_CUTOFF", param.TRIOZON_A+5*param.TRIOZON_DELTA)
param.createdefault("TRIOZON_Z_CUTOFF", (param.TRIOZON_CUTOFF**2 - (0.95*param.GRAPHITE_INTERLAYER_DISTANCE)**2)**0.5)
