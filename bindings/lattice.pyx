# distutils: language = c++

from Lattice cimport (
    Lattice,                                # Base Class
    Hexagonal, Square, Triangular, UJack,   # Regular2d
    Cubic, BCC,                             # Regular3d
    Cubic4d,                                # Regular4d
    FixedZ,                                 # Random coordinated
    SquRand, TriRand, CubRand,              # Random bonds
    SWNetwork,                              # Small world
    Point2d, Point3d, Point4d,              # Helpers
)


# cname hack for non-type templates.
# see https://stackoverflow.com/questions/53582945/wrapping-c-code-with-function-pointer-as-template-parameter-in-cython
cdef extern from *:
    ctypedef size_t T4 "4"
    ctypedef size_t T3 "3"
    ctypedef size_t T2 "2"
    ctypedef size_t T1 "1"
    ctypedef size_t T0 "0"


cdef class __LatticeMixin:
    """ All the python methods that come with being a Lattice. 
    
    This cannot be used directly. Classes derived from Lattice that are to be
    exposed to Python must derive from this class and provide a __cinit__ method
    that creates and stores an instance of the lattice type. 

    See TriangularLattice for an example.
    
    """

    cdef Lattice *_p

    def __dealloc__(self):
        del self._p

    def getLength(self):
        return self._p.getLength()

    def getNumSites(self):
        return self._p.getNumSites()

    def getNumActive(self):
        return self._p.getNumActive()

    def getNumPresent(self):
        return self._p.getNumPresent()

    def getNumClusters(self):
        return self._p.getNumClusters()

    def getMaxMass(self):
        return self._p.getMaxMass()

    def isActive(self, i):
        return self._p.isActive(i)

    def isPresent(self, i):
        return self._p.isPresent(i)

    def setActiveLevel(self, i, level):
        return self._p.setActiveLevel(i, level)

    def setPresentLevel(self, i, level):
        return self._p.setPresentLevel(i, level)

    def getClusterLabel(self, i):
        return self._p.getClusterLabel(i)

    def getClusterSize(self, c):
        return self._p.getClusterSize(c)

    def isSpanning(self, c):
        return self._p.isSpanning(c)

    def getNumNeighbors(self, i):
        return self._p.getNumNeighbors(i)

    def getNbr(self, i, j):
        return self._p.getNbr(i, j)

    def getNumActiveNeighbors(self, i):
        return self._p.getNumActiveNeighbors(i)

    def getDims(self):
        return self._p.getDims()

    # These are to make the Python code easier to use

    def getCoords(self, i):
        """ Get the coordinates coresponding to a lattice point. 

        Any lattice that is not 2, 3, or 4 dimensions is treated as
        2-dimeinsional.
        
        """
        dims = self._p.getDims()
        if dims == 4:
            p4 = new Point4d(i, self._p.getLength())
            t = (p4.w, p4.x, p4.y, p4.z)
            del p4
        elif dims == 3:
            p3 = new Point3d(i, self._p.getLength())
            t = (p3.x, p3.y, p3.z)
            del p3
        else:
            p2 = new Point2d(i, self._p.getLength())
            t = (p2.x, p2.y)
        return t

    def iter_nbrs(self, i):
        n = self._p.getNumNeighbors(i)
        for j in range(n):
            yield self._p.getNbr(i, j)

    def __len__(self):
        return self._p.getNumSites()

    def __iter__(self):
        for i in range(self._p.getNumSites()):
            yield i


# 2D Regular Lattices

cdef class HexagonalLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new Hexagonal[T0](len_, 0)
        elif pbc == 1:
            self._p = new Hexagonal[T1](len_, 0)
        else:
            self._p = new Hexagonal[T2](len_, 0)


cdef class SquareLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new Square[T0](len_, 0)
        elif pbc == 1:
            self._p = new Square[T1](len_, 0)
        else:
            self._p = new Square[T2](len_, 0)


cdef class TriangularLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new Triangular[T0](len_, 0)
        elif pbc == 1:
            self._p = new Triangular[T1](len_, 0)
        else:
            self._p = new Triangular[T2](len_, 0)


cdef class UJackLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new UJack[T0](len_, 0)
        elif pbc == 1:
            self._p = new UJack[T1](len_, 0)
        else:
            self._p = new UJack[T2](len_, 0)


cdef class CubicLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new Cubic[T0](len_, 0)
        elif pbc == 1:
            self._p = new Cubic[T1](len_, 0)
        elif pbc == 2:
            self._p = new Cubic[T2](len_, 0)
        else:
            self._p = new Cubic[T3](len_, 0)


cdef class BCCLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new BCC[T0](len_, 0)
        elif pbc == 1:
            self._p = new BCC[T1](len_, 0)
        elif pbc == 2:
            self._p = new BCC[T2](len_, 0)
        else:
            self._p = new BCC[T3](len_, 0)


cdef class Cubic4dLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, unsigned short pbc):
        if pbc == 0:
            self._p = new Cubic4d[T0](len_, 0)
        elif pbc == 1:
            self._p = new Cubic4d[T1](len_, 0)
        elif pbc == 2:
            self._p = new Cubic4d[T2](len_, 0)
        elif pbc == 3:
            self._p = new Cubic4d[T3](len_, 0)
        else:
            self._p = new Cubic4d[T4](len_, 0)


cdef class FixedZLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, size_t z):
        self._p = new FixedZ(len_, 0, z)


cdef class SquRandLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, double fb):
        self._p = new SquRand(len_, 0, fb)


cdef class TriRandLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, double fb):
        self._p = new TriRand(len_, 0, fb)


cdef class CubRandLattice(__LatticeMixin):

    def __cinit__(self, size_t len_, double fb):
        self._p = new CubRand(len_, 0, fb)


cdef class SWNLattice(__LatticeMixin):

    def __cinit__(self, size_t nsites, size_t dim, double alpha, double gamma):
        self._p = new SWNetwork(nsites, 0, dim, alpha, gamma)