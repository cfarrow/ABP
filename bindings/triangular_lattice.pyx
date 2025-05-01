# distutils: language = c++

from Triangular cimport Triangular

cdef class TriangularLattice:
    cdef Triangular* _p;

    def __cinit__(self, size_t len_, size_t id_=0, short pbc=2):
        self._p = new Triangular(len_, id_, pbc)

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

    def iter_nbrs(self, i):
        n = self._p.getNumNeighbors(i)
        for j in range(n):
            yield self._p.getNbr(i, j)