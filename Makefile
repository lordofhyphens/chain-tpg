# CUDD compilation
CUDD:=./cudd
CUDDLIBS=obj cudd epd mtr util st 
CUDDLIB_SEARCH=$(foreach libs,$(CUDDLIBS),-L$(CUDD)/$(libs)) 
CUDDLIB_FLAGS=$(foreach libs,$(CUDDLIBS),-l$(libs))
CUDDFLAGS:=-I$(CUDD)/include

# CPPUTEST 
CPPUTEST_FLAGS:=-I$(CPPUTEST_HOME)/include
CPPUTEST_LIBS:=-lCppUTest -lCppUTestExt

objs=getpis.o bdd_img.o explore.o cudd_ckt.o
TEST=$(foreach test,$(objs:.o=_test.cpp) AllTests.cpp,tests/${test})
CPPFLAGS:=$(CUDDFLAGS) -Iutil -D_GLIBCXX_PARALLEL 
CXXFLAGS:= -O3 -mtune=native -DHAVE_IEEE_754 -DBSD -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 -DCPU -g -Wall -std=c++11 -march=native -fopenmp 
LDFLAGS+=-static $(CUDDLIB_SEARCH) -L./util 

LIBS:=$(CUDDLIB_FLAGS) -lcktutil -lm -lz
CXX=g++
.PHONY: all test clean
.SUFFIXES: cc c cpp

explore: $(objs)
	echo $(LIBS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

explore.o: explore.cpp cudd_ckt.h bdd_img.h 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
	
bdd_img.o: bdd_img.cpp bdd_img.h cudd_ckt.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@ $(LIBS)
cudd_ckt.o: cudd_ckt.cpp cudd_ckt.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

util/libcktutil.a: util/ckt.cc
	git submodule init
	git submodule update
	cd util && $(MAKE) libs

$(CUDD)/cudd/libcudd.a:	$(CUDD)/Makefile
	cd cudd && $(MAKE)

$(CUDD)/Makefile:
	git submodule init
	git submodule update
clean: 
	rm -f explore *.o tests/*.o testsuite

AllTests.o: $(TEST)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $^ -o $@

testsuite: $(subst explore.o,,$(objs)) $(TEST:.cpp=.o)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CPPUTEST_FLAGS) $(LDFLAGS) -L${CPPUTEST_HOME}/lib $^ $(LIBS) $(CPPUTEST_LIBS) -o $@ 

test: testsuite
	./$<
