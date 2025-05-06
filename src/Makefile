# CXXFLAGS
OPTFLAGS = -fPIC -O3 -ffast-math -funroll-loops #-fno-inline 
DBGFLAGS = -g -Wall -Wextra

ifdef OPT
CXXFLAGS = $(OPTFLAGS)
else
CXXFLAGS = $(DBGFLAGS)
endif


# Rules of library files
OBJS = ABPprocs.o BCC.o Cubic4d.o Cubic.o CubRand.o FixedZ.o Hexagonal.o Square.o SquRand.o SWNetwork.o Triangular.o TriRand.o UJack.o Lattice.o LatticeRegular.o LatticeRegular2d.o LatticeRegular3d.o LatticeRegular4d.o LatticeRandom.o
# 2D lattices
D2OBJS = Hexagonal.o Square.o Triangular.o UJack.o
# 3D lattices
D3OBJS = BCC.o Cubic.o
# 4D lattices
D4OBJS = Cubic4d.o
# Random Lattices
RANDOBJS = CubRand.o FixedZ.o SquRand.o SWNetwork.o TriRand.o

# program rules
%: %.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $*.out

# library rules
$(OBJS): MersenneTwister.h Point.h Lattice.h Lattice.cpp
$(D2OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular2d.h LatticeRegular2d.cpp
$(D3OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular3d.h LatticeRegular3d.cpp
$(D4OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular4d.h LatticeRegular4d.cpp
$(RANDOBJS): LatticeRandom.h LatticeRandom.cpp
# All other libraries
%.o: %.cpp %.h
	$(CXX) $< $(CXXFLAGS) -c -o $@

# Utility
.PHONY : clean
clean :
	-rm *.o

# Include Make rules for tests. Edit tests.mk when adding tests.
include tests.mk
