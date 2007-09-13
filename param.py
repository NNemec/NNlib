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
    def __init__(self):
        self.__dict__['defaults'] = {}
        self.__dict__['temp'] = {}
        self.__dict__['data'] = DataDummy()

    def setdefaults(self):
        for name in self.data._v_attrnamesuser:
            delattr(self.data,name)

        for name,value in self.defaults.iteritems():
            setattr(self.data,name,value)

    def createdefault(self,name,value):
        self.defaults[name] = value
        if not hasattr(self.data,name):
            setattr(self.data,name,value)

    def write(self,data=None):
	import tables as pytables

        if isinstance(data,pytables.Group):
            data = data._v_attrs
        elif isinstance(data,pytables.Leaf):
            data = data.attrs
        for name in data._v_attrnamesuser:
            delattr(data,name)
        for name in self.data._v_attrnamesuser:
            setattr(data,name,getattr(self.data,name))
        self.__dict__['data'] = data

    def reset(self):
	self.write(DataDummy())

    def read(self,data):
	import tables as pytables

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


# param.createdefault("LATTICE_CONSTANT", param.GRAPHENE_CC_DISTANCE)
# param.createdefault("HOPPING_CONSTANT", param.GRAPHENE_1STNN_HOPPING)


