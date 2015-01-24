
# CUDD compilation stuff
CUDD:=./cudd
CUDDLIBS=obj cudd epd mtr util st 
CUDDLIB_SEARCH=$(foreach libs,$(CUDDLIBS),-L$(CUDD)/$(libs)) 
CUDDLIB_FLAGS=$(foreach libs,$(CUDDLIBS),-l$(libs))


objs=bdd_img.o explore.o cudd_ckt.o
CUDDFLAGS:=-I$(CUDD)/include
CPPUTEST_FLAGS:=-I$(CPPUTEST_HOME)/include

CXXFLAGS:=$(CUDDFLAGS) $(CPPUTEST_FLAGS) -O3 -mtune=native -DHAVE_IEEE_754 -DBSD -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 -DCPU -g -Wall -std=c++11 -march=native -fopenmp -D_GLIBCXX_PARALLEL 
LDFLAGS+=$(CUDDLIB_SEARCH) -L${CPPUTEST_HOME}/lib -L./util 

LIBS:=$(CUDDLIB_FLAGS) -lcktutil -lm -lz -lCppUTestExt -lCppUTest
CXX=g++
.PHONY: all test

explore: $(objs)
	echo $(LIBS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

explore.o: explore.cc cudd_ckt.h bdd_img.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
bdd_img.o: bdd_img.cc bdd_img.h cudd_ckt.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@ $(LIBS)
cudd_ckt.o: cudd_ckt.cc cudd_ckt.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

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
	rm -f explore *.o
AllTests.o: AllTests.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

testsuite: AllTests.o bdd_img.o getpis_test.o 
	$(CXX)  $(LDFLAGS) $(CXXFLAGS) $^ $(LIBS) -o $@ 
test: testsuite
	./$@
