# distutils: language = c++


cdef extern from "Triangular.cpp":
    pass

cdef extern from "LatticeRegular2d.cpp":
    pass

cdef extern from "LatticeRegular.cpp":
    pass

cdef extern from "Lattice.cpp":
    pass


cdef extern from "Triangular.h":
    cdef cppclass Triangular:
        Triangular(size_t len_, size_t id_, short pbc) except +
        size_t getLength()
        size_t getNumSites()
        size_t getNumActive()
        size_t getNumPresent()
        size_t getNumClusters()
        size_t getMaxMass()

        #/* Functions that act on single sites */
        bint isActive(size_t)
        bint isPresent(size_t)
        void setActiveLevel(size_t, bint)
        void setPresentLevel(size_t, bint)
        size_t getClusterLabel(size_t i)
        size_t getClusterSize(size_t i)
        bint isSpanning(size_t)
        size_t getNumNeighbors(size_t) 
        size_t getNbr(size_t, size_t)
        size_t getNumActiveNeighbors(size_t)
