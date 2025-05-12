# distutils: language = c++


cdef extern from "rand.cpp":
    pass

cdef extern from "Lattice.cpp":
    pass

cdef extern from "LatticeRegular.cpp":
    pass

cdef extern from "LatticeRegular2d.cpp":
    pass

cdef extern from "Hexagonal.cpp":
    pass

cdef extern from "Square.cpp":
    pass

cdef extern from "Triangular.cpp":
    pass

cdef extern from "UJack.cpp":
    pass

cdef extern from "LatticeRegular3d.cpp":
    pass
    
cdef extern from "Cubic.cpp":
    pass
    
cdef extern from "BCC.cpp":
    pass

cdef extern from "LatticeRegular4d.cpp":
    pass

cdef extern from "Cubic4d.cpp":
    pass

cdef extern from "LatticeRandom.cpp":
    pass

cdef extern from "FixedZ.cpp":
    pass

cdef extern from "FixedZ.cpp":
    pass

cdef extern from "SquRand.cpp":
    pass

cdef extern from "TriRand.cpp":
    pass

cdef extern from "CubRand.cpp":
    pass

cdef extern from "SWNetwork.cpp":
    pass


cdef extern from "Point.h":
    cdef cppclass Point2d:
        Point2d()
        void update(size_t index, size_t len_) nogil
        size_t x, y

    cdef cppclass Point3d:
        Point3d()
        void update(size_t index, size_t len_) nogil
        size_t x, y, z

    cdef cppclass Point4d:
        Point4d()
        void update(size_t index, size_t len_) nogil
        size_t w, x, y, z


cdef extern from "Lattice.h":
    cdef cppclass Lattice:
        Lattice(size_t len_, size_t id_)
        
        void activateSites() nogil 
        void labelClusters(size_t) nogil
        size_t getDims() nogil
        size_t getLength() nogil
        size_t getNumSites() nogil
        size_t getNumActive() nogil
        size_t getNumPresent() nogil
        size_t getNumClusters() nogil
        size_t getMaxMass() nogil

        #/* Functions that act on single sites */
        bint isActive(size_t) nogil
        bint isPresent(size_t) nogil
        void setActiveLevel(size_t, bint) nogil
        void setPresentLevel(size_t, bint) nogil
        size_t getClusterLabel(size_t i) nogil
        size_t getClusterSize(size_t i) nogil
        bint isSpanning(size_t) nogil
        size_t getNumNeighbors(size_t)  nogil
        size_t getNbr(size_t, size_t) nogil
        size_t getNumActiveNeighbors(size_t) nogil


cdef extern from "Hexagonal.h":
    cdef cppclass Hexagonal[T](Lattice):
        Hexagonal(size_t len_, size_t id_)

cdef extern from "Square.h":
    cdef cppclass Square[T](Lattice):
        Square(size_t len_, size_t id_)

cdef extern from "Triangular.h":
    cdef cppclass Triangular[T](Lattice):
        Triangular(size_t len_, size_t id_)

cdef extern from "UJack.h":
    cdef cppclass UJack[T](Lattice):
        UJack(size_t len_, size_t id_)

cdef extern from "Cubic.h":
    cdef cppclass Cubic[T](Lattice):
        Cubic(size_t len_, size_t id_)

cdef extern from "BCC.h":
    cdef cppclass BCC[T](Lattice):
        BCC(size_t len_, size_t id_)

cdef extern from "Cubic4d.h":
    cdef cppclass Cubic4d[T](Lattice):
        Cubic4d(size_t len_, size_t id_)

cdef extern from "FixedZ.h":
    cdef cppclass FixedZ(Lattice):
        FixedZ(size_t len_, size_t id_, size_t z)

cdef extern from "SquRand.h":
    cdef cppclass SquRand(Lattice):
        SquRand(size_t len_, size_t id_, double fb)

cdef extern from "TriRand.h":
    cdef cppclass TriRand(Lattice):
        TriRand(size_t len_, size_t id_, double fb)

cdef extern from "CubRand.h":
    cdef cppclass CubRand(Lattice):
        CubRand(size_t len_, size_t id_, double fb)

cdef extern from "SWNetwork.h":
    cdef cppclass SWNetwork(Lattice):
        SWNetwork(size_t nsites, size_t id_, size_t dim, double alpha, double gamma)
