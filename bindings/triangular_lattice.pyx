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