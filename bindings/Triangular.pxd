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
        #size_t getNumSites() {return num_sites;}
        #size_t getNumActive() {return num_active_sites;}
        #size_t getNumPresent() {return num_present_sites;}
        #size_t getNumClusters() {return Clusters.size();}
        #size_t getMaxMass() {return (Clusters.empty() ? 0 : *(max_element(Clusters.begin(), Clusters.end()))); }

        #/* Functions that act on single sites */
        #bool isActive(size_t);
        #bool isPresent(size_t);
        #/* the following do not do any checking. It is up to the algorithm. 'if' statements slow things down! */
        #void setActiveLevel(size_t, bool);
        #void setPresentLevel(size_t, bool);
        #/* get cluster label of site */
        #size_t getClusterLabel(size_t i) {return cluster_label[i]; }
        #/* get size of cluster, as referenced by its label */
        #size_t getClusterSize(size_t i) {return Clusters[i]; }

        #/* virtual functions */
        #virtual bool isSpanning(size_t = 1);

        #/* virtual function - see notes below */
        #virtual size_t getNumNeighbors(size_t); 
        #virtual size_t getNbr(size_t, size_t);
        #virtual size_t getNumActiveNeighbors(size_t);
