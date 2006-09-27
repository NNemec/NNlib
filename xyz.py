from calc import *
from param import param
from copy import deepcopy
#from cnt import armchair, zigzag, chiral, swcnt

class atom:
    def __init__(self,typ,pos,rot=None):
        self.typ = typ
        self.pos = asarray(pos)
        assert self.pos.shape == (3,)
	if rot is not None:
	    self.rot = asmatrix(rot)
	    assert self.rot.shape == (3,3)

    def shift(self,disp):
        return atom(self.typ,self.pos+disp)

    def rotate(self,rot):
	rot = asmatrix(rot)
        assert rot.shape == (3,3)
	pos = dot(rot,self.pos)
	if hasattr(self,rot):
	    rot = dot(rot,self.rot)
	else:
	    rot = None
        return atom(self.typ,pos=pos,rot=rot)

class structure:
    def __init__(self):
        self.atoms = []

    def shift(self,disp):
        res = deepcopy(self)
        for i in range(len(res.atoms)):
            res.atoms[i] = res.atoms[i].shift(disp)
        return res

    def rotate(self,rot):
        res = deepcopy(self)
        for i in range(len(res.atoms)):
            res.atoms[i] = res.atoms[i].rotate(rot)
        return res

    def sorted_by_z(self):
        res = deepcopy(self)
        res.atoms = sorted(self.atoms,key=lambda a:a.pos[2])
        return res

    def write_xyz_file(self,filename):
        f = file(filename,'w')
        print >> f, len(self.atoms), " AN"
        print >> f, "xyz-file created by Norbert Nemec (Python)"
        for at in self.atoms:
            print >> f, at.typ, "\t", at.pos[0], "\t", at.pos[1], "\t", at.pos[2]
        f.close()

    def read_xyz_file(self,filename):
        f = file(filename,'r')
        Natoms = int(f.readline().strip().split()[0])
        f.readline()
        for n in range(Natoms):
            typ,x,y,z = f.readline().strip().split()
            self.atoms.append(atom(typ,[float(x),float(y),float(z)]))
        assert f.read().strip() == ''
        f.close()

class chain(structure):
    def __init__(self,period):
        structure.__init__(self)
        self.period = asarray(period)
        assert self.period.shape == (3,)

    def multiply(self,N):
        res = chain(N*self.period)
        for n in range(N):
            for atom in self.atoms:
                res.atoms.append(atom.shift(n*self.period))
        return res

class sheet(structure):
    def __init__(self,period):
        structure.__init__(self)
        assert len(period) == 2
        a = [asarray(period[0]),asarray(period[1])]
        assert a[0].shape == (3,)
        assert a[1].shape == (3,)
        self.period = a
        ez = [0,0,1]
        spade = dot(a[0],cross(a[1],ez))
        assert spade != 0.
        self.rzp = [
            cross(a[1],ez) / spade,
            cross(ez,a[0]) / spade,
        ]

    def multiply(self,N0,N1):
        res = sheet([N0*self.period[0],N1*self.period[1]])
        for n0 in range(N0):
            for n1 in range(N1):
                for atom in self.atoms:
                    res.atoms.append(atom.shift(n0*self.period[0] + n1*self.period[1]))
        return res

def square_ladder(N):
    spacing = param.GRAPHENE_CC_DISTANCE
    res = chain((0,0,spacing))
    for n in range(N):
        res.atoms.append(atom('C',(spacing*(n-(N-1)*0.5),0,0)),rot=eye(3))
    return res

def square_tube(N):
    spacing = param.GRAPHENE_CC_DISTANCE
    radius = N*spacing/(2*pi)
    res = chain((0,0,spacing))
    for n in range(N):
        phi = 2*pi*n/N
        res.atoms.append(atom('C',(cos(phi)*radius,sin(phi)*radius,0)),rot=[[cos(phi),-sin(phi),0],[sin(phi),cos(phi),0],[0,0,1]])
    return res

def linchain():
    return square_ladder(1)

def merge(xyz_A,xyz_B):
    if isinstance(xyz_A,chain):
        assert isinstance(xyz_B,chain)
        assert((xyz_A.period == xyz_B.period).all())
        res = chain(xyz_A.period)
        res.atoms = xyz_A.atoms + xyz_B.atoms
    elif isinstance(xyz_A,sheet):
        assert isinstance(xyz_B,sheet)
        assert((xyz_A.period[0] == xyz_B.period[0]).all())
        assert((xyz_A.period[1] == xyz_B.period[1]).all())
        res = sheet(xyz_A.period)
        res.atoms = xyz_A.atoms + xyz_B.atoms
    return res

class aperiodic:
    def __init__(self,length,cellsize):
        self.length = length
        self.cellsize = cellsize
        self.cellcount = int(-((-length)//cellsize))
        self.cells = []
        for i in range(self.cellcount):
            self.cells += [structure()]

    def add_chain(self,chain):
        ch = chain.sorted_by_z()
        curr_cell = -1
        curr_cell_boundary = 0.0
        while curr_cell < self.cellcount:
            for at in ch.atoms:
                if at.pos[2] >= curr_cell_boundary:
                    curr_cell += 1
                    curr_cell_boundary += self.cellsize
                if curr_cell >= 0 and curr_cell < self.cellcount:
                    self.cells[curr_cell].atoms.extend([at])
            ch = ch.shift(ch.period)

    def add_structure(self,structure):
        ch = structure.sorted_by_z()
        curr_cell = -1
        curr_cell_boundary = 0.0
        for at in ch.atoms:
            if at.pos[2] >= curr_cell_boundary:
                curr_cell += 1
                curr_cell_boundary += self.cellsize
            if curr_cell >= 0 and curr_cell < self.cellcount:
                self.cells[curr_cell].atoms.extend([at])

    def write_xyz_file(self,filename,charge=None):
        f = file(filename,'w')
        print >> f, sum([len(c.atoms) for c in self.cells]), " AN"
        print >> f, "xyz-file created by Norbert Nemec (Python)"
        if charge is None:
            for c in self.cells:
                for at in c.atoms:
                    print >> f, at.typ + "\t%f\t%f\t%f"%(tuple(at.pos))
        else:
            for c,ch in zip(self.cells,charge):
                for at,chg in zip(c.atoms,ch):
                    print >> f, at.typ + "\t%f\t%f\t%f"%(tuple(at.pos)) + "\t%g"%chg
        f.close()

if __name__ == "__main__":
    param.setdefaults
    dummy = square_ladder(3)
    dummy = dummy.multiply(3)
    for atom in dummy.atoms:
        print atom.typ, ': ', atom.pos
