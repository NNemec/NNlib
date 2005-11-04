from calc import *
from param import param
from copy import deepcopy

class atom:
    def __init__(self,typ,pos):
        self.typ = typ
        assert pos.shape == (3,)
        self.pos = pos

    def shift(self,disp):
        assert disp.shape == (3,)
        return atom(self.typ,self.pos+disp)

    def rotate(self,rot):
        assert rot.shape == (3,3)
        return atom(self.typ,dot(rot,self.pos))

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
            self.atoms.append(atom(typ,c_[float(x),float(y),float(z)]))
        assert f.read().strip() == ''
        f.close()

class chain(structure):
    def __init__(self,period):
        structure.__init__(self)
        self.period = period

    def multiply(self,N):
        res = chain(N*self.period)
        for n in range(N):
            for atom in self.atoms:
                res.atoms.append(atom.shift(n*self.period))
        return res

def square_ladder(N,spacing=1):
    spacing = param.GRAPHENE_CC_DISTANCE
    res = chain(c_[0,0,spacing])
    at = atom('C',c_[0,0,0])
    sh = c_[spacing,0,0]
    for n in range(N):
        res.atoms.append(at.shift(sh*n))
    return res

def linchain():
    return square_ladder(1)

def armchair(N):
    CC_distance = param.GRAPHENE_CC_DISTANCE

    period = c_[0,0,sqrt(3)*CC_distance]
    circumference = N*2*1.5*CC_distance
    radius = circumference/(2*pi)
    angle = 2*pi/(6*N)

    rot = Matrix(eye(3))
    rot[0,0] = cos(angle)
    rot[1,1] = cos(angle)
    rot[0,1] = sin(angle)
    rot[1,0] = -sin(angle)
    rot2 = rot*rot
    rot3 = rot2*rot
    rot4 = rot2*rot2
    rot6 = rot3*rot3

    at = atom('C',c_[radius,0,0])

    cell = structure()
    cell.atoms = [
        at.shift(period/2),
        at.rotate(rot),
        at.rotate(rot3),
        at.rotate(rot4).shift(period/2)
    ]

    res = chain(period)
    for i in range(N):
        res.atoms.extend(cell.atoms)
        cell = cell.rotate(rot6)
    res.radius = radius

    return res

def zigzag(N):
    CC_distance = param.GRAPHENE_CC_DISTANCE

    period = c_[0,0,3*CC_distance]
    circumference = N*sqrt(3)*CC_distance
    radius = circumference/(2*pi)
    angle = pi/N

    rot = Matrix(eye(3))
    rot[0,0] = cos(angle)
    rot[1,1] = cos(angle)
    rot[0,1] = sin(angle)
    rot[1,0] = -sin(angle)
    rot2 = rot*rot

    at = atom('C',c_[radius,0,0])

    cell = structure()
    cell.atoms = [
        at.shift(period/6),
        at.shift(period/2),
        at.rotate(rot),
        at.rotate(rot).shift(period/1.5)
    ]

    res = chain(period)
    for i in range(N):
        res.atoms.extend(cell.atoms)
        cell = cell.rotate(rot2)
    res.radius = radius

    return res

def swcnt(V):
    assert len(V)==2
    if V[0] == V[1]:
        return armchair(V[0])
    elif V[1] == 0:
        return zigzag(V[0])
    elif V[0] == 0:
        return zigzag(V[1])
    else:
        raise RuntimeError, "chiral tubes not implemented"

def merge(chain_A,chain_B):
    assert(chain_A.period == chain_B.period)
    res = chain(chain_A.period)
    res.atoms = chain_A.atoms + chain_B.atoms
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
