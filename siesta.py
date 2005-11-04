from calc import *
from units import *

class readCUBE:
    def __init__(self,filename):
	def cvlen(s):
	    v = float(s)
	    if v < 0:
                return -v*bohr
            else:
                return v*angstrom
		
        f = open(filename)
        f.readline()
        f.readline()
        ll = f.readline().split()
        self.Natom = int(ll[0])
        self.origin = array((cvlen(ll[1]),cvlen(ll[2]),cvlen(ll[3])))
        self.voxels = [None]*3
        self.vector = [None]*3
        for i in range(3):
            ll = f.readline().split()
            self.voxels[i] = int(ll[0])
            self.vector[i] = array((cvlen(ll[1]),cvlen(ll[2]),cvlen(ll[3])))
        self.atomType = [None]*self.Natom
        self.atomPos = [None]*self.Natom
        for n in range(self.Natom):
            ll = f.readline().split()
            self.atomType[n] = int(ll[0])
	    dummy = float(ll[1])
            self.atomPos[n] = array((cvlen(ll[2]),cvlen(ll[3]),cvlen(ll[4])))
        data = f.read().split()
        assert len(data) == prod(self.voxels)
        self.data = zeros(self.voxels,'d')
        i = 0
        for ix in range(self.voxels[0]):
         for iy in range(self.voxels[1]):
          for iz in range(self.voxels[2]):
            self.data[ix,iy,iz] = float(data[i])
            i += 1

    def writeCUBE(self,filename):
        f = open(filename,'w')
        print >>f," comment"
        print >>f," comment"
	print >>f,"% 5i% 12.6f% 12.6f% 12.6f"%(self.Natom,self.origin[0],self.origin[1],self.origin[2])

        for i in range(3):
    	    print >>f,"% 5i% 12.6f% 12.6f% 12.6f"%(self.voxels[i],self.vector[i][0],self.vector[i][1],self.vector[i][2])

        for n in range(self.Natom):
    	    print >>f,"% 5i% 12.6f% 12.6f% 12.6f% 12.6f"%(
		self.atomType[n],0.0,
		self.atomPos[n][0],
		self.atomPos[n][1],
		self.atomPos[n][2],
	    )

	l = 0
	for d in ravel(self.data):
	    f.write('  %.5E'%d)
	    l += 1
	    if l == 6:
		f.write('\n')
		l = 0
	if l > 0:
	    f.write('\n')
	f.close()	

    def writeVTK(self,filename,name='cubedata'):
        import pyvtk

        point_data = pyvtk.PointData(
	    pyvtk.Scalars(
		reshape(transpose(self.data),(-1,)),
		name=name
	    )
	)

        assert not self.vector[0][1]
        assert not self.vector[0][2]
        assert not self.vector[1][2]
        assert not self.vector[1][0]
        assert not self.vector[2][0]
        assert not self.vector[2][1]

        grid = pyvtk.StructuredPoints(
	    self.voxels, 
	    self.origin,
	    (
		self.vector[0][0],
		self.vector[1][1],
        	self.vector[2][2],
    	    ),
	)
        vtkdata = pyvtk.VtkData(grid, point_data)
        vtkdata.tofile(filename,format='binary')

    def writeNetCDF(self,filename):
	from Scientific.IO.NetCDF import NetCDFFile
	ncfile = NetCDFFile(filename, 'w')
	
	for dim,i in (('x',0),('y',1),('z',2)):
	    ncfile.createDimension(dim,self.voxels[i])
	    ncfile.createVariable(dim,'d',(dim,))[:] = arange(self.voxels[i])*self.vector[i][i]
	
	ncfile.createVariable('data','d',('x','y','z'))
	ncfile.variables['data'][:] = self.data

	ncfile.Natom = self.Natom
	ncfile.origin = self.origin
	for n in range(self.Natom):
	    setattr(ncfile,'atom%i.Type'%n,self.atomType[n])
	    setattr(ncfile,'atom%i.Pos'%n,self.atomPos[n])

	ncfile.close()


class readBANDS:
    def __init__(self,filename):
        f = file(filename)

        ll = f.readline().split()
        self.E_F = float(ll[0]) * eV

        ll = f.readline().split()
        self.k_min,self.k_max = float(ll[0]),float(ll[1])

        ll = f.readline().split()
        self.E_min,self.E_max = float(ll[0]),float(ll[1]) * eV

        ll = f.readline().split()
        self.Nbands,self.Nspin,self.Nk = int(ll[0]),int(ll[1]),int(ll[2])

        ll = f.read().split()

        self.x = zeros((self.Nk),'d')
        self.E = zeros((self.Nk,self.Nbands),'d')

        for k in range(self.Nk):
            self.x[k] = float(ll[k*(self.Nbands+1)])
            for b in range(self.Nbands):
                self.E[k,b] = float(ll[k*(self.Nbands+1)+b+1]) * eV

        ll = ll[self.Nk*(self.Nbands+1):]
        self.Npoint = int(ll[0])
        ll = ll[1:]

        self.pointX = [None]*self.Npoint
        self.pointName = [None]*self.Npoint
        for p in range(self.Npoint):
            self.pointX[p] = float(ll[2*p])
            self.pointName[p] = ll[2*p+1][1:-1]

class readDOS:
    def __init__(self,filename):
        f = open(filename)

        f.readline()
        f.readline()
        f.readline()
        f.readline()

        ll = f.readline().split()
        self.Nbands = int(ll[-3])
        self.Nspin = int(ll[-2])
        self.Nk = int(ll[-1])

        ll = f.readline().split()
        self.E_F = float(ll[-2]) * eV

        ll = f.readline().split()
        self.broadening = float(ll[-2]) * eV

        ll = f.readline().split()
        self.N_e1 = float(ll[-4]) * eV
        self.N_e2 = float(ll[-3]) * eV

        f.readline()
        f.readline()
        f.readline()

        lines = f.readlines()

        self.N = zeros(len(lines),'d')
        self.E = zeros(len(lines),'d')

        for i in range(len(lines)):
            ll = lines[i].split()
            self.N[i] = float(ll[0])
            self.E[i] = float(ll[1]) * eV

def eig2dos(basename,broadening,steps,E_min,E_max):
    from glob import glob
    eigname = glob(basename+'.d/*.EIG')
    assert len(eigname) == 1
    eigname, = eigname
    dosname = basename+'.DOS'
    fin = open(eigname)
    fout = open(dosname,'w')

    from subprocess import Popen, PIPE
    p = Popen('eig2dos', shell=True, stdin=PIPE, stdout=fout)

    p.stdin.write(fin.readline())
    p.stdin.write("%f %i %f %f\n"%(broadening/u.eV,steps,E_min/u.eV,E_max/u.eV))
    d = fin.read()
    p.stdin.write(d)
    p.stdin.close()
    p.wait()
    fin.close()

def rho2cube(basename):
    from glob import glob
    SystemLabel = glob(basename+'.d/*.RHO')[0].split('/')[-1][:-4]
    from subprocess import Popen, PIPE
    p = Popen('grid2cube', shell=True, stdin=PIPE, cwd=basename+'.d')

    p.stdin.write("""\
%(SystemLabel)s
rho
0.0 0.0 10.0
1
unformatted
"""%locals())
    p.stdin.close()
    p.wait()

def ldos2cube(basename):
    from glob import glob
    SystemLabel = glob(basename+'.d/*.RHO')[0].split('/')[-1][:-4]
    from subprocess import Popen, PIPE
    p = Popen('grid2cube', shell=True, stdin=PIPE, cwd=basename+'.d')

    p.stdin.write("""\
%(SystemLabel)s
ldos
0.0 0.0 10.0
1
unformatted
"""%locals())
    p.stdin.close()
    p.wait()

