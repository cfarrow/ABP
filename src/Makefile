# CFLAGS
OPTFLAGS = -fPIC -O3 -ffast-math -funroll-loops #-fno-inline 
#OPTFLAGS = -O3 -funroll-loops #-fno-inline 
DBGFLAGS = -gstabs+ -Wall

CC = g++
LINK = g++

ifdef DEBUG
CFLAGS = $(DBGFLAGS)
else
CFLAGS = $(OPTFLAGS)
endif


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
	$(LINK) $(CFLAGS) $^ -o $*

# library rules
$(OBJS): MersenneTwister.h Lattice.h Lattice.cpp
$(D2OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular2d.h LatticeRegular2d.cpp
$(D3OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular3d.h LatticeRegular3d.cpp
$(D4OBJS): LatticeRegular.h LatticeRegular.cpp LatticeRegular4d.h LatticeRegular4d.cpp
$(RANDOBJS): LatticeRandom.h LatticeRandom.cpp
# All other libraries
%.o: %.cpp %.h
	$(CC) $< $(CFLAGS) -c -o $@

# Utility
.PHONY : clean
clean :
	-rm *.o


# $*
#     The stem with which an implicit rule matches. If the target is
#     'dir/a.foo.b' and the target pattern is 'a.%.b' then the stem is
#     'dir/foo'. The stem is useful for constructing names of related files.
# 
# $@
#     The file name of the target of the rule.
# 
# $<
#     The name of the first dependency.
# 
# $?
#     The names of all the dependencies that are newer than the target, with
#     spaces between them.
# 
# $^
#     The names of all the dependencies, with spaces between them.
