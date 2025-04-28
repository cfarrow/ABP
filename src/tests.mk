# This is not meant to be run on it's own. It is included in the main Makefile.
#   make [tests]  - makes everything.
#   make TARGET - makes the given target.

# All tests produced by this Makefile.  Add new tests here and deps below.
TESTS = test_wordops

test_wordops.o : MersenneTwister.h Lattice.h Lattice.cpp

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = /usr/src/googletest/googletest

TEST_DIR = ./tests

# Extra flags passed to the preprocessor for testing.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
TESTFLAGS = -isystem $(GTEST_DIR)/include


# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

tests : $(TESTS)

.PHONY : clean_tests
clean_tests :
	rm -f $(TESTS) gtest.a gtest_main.a

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(TESTFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) $(TESTFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(TESTFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) $(TESTFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^


# Build any test. Add exceptions as needed.

test_%.o : $(TEST_DIR)/test_%.cpp $(GTEST_HEADERS)
	$(CXX) $< $(CXXFLAGS) $(TESTFLAGS) -c -o $@

test_% : test_%.o $(OBJS) gtest_main.a
	$(CXX) $(CXXFLAGS) -lpthread $^ -o $@
