CUDD:=./cudd
CUDDLIBS=cudd mtr util dddmp epd obj st 
objs=bdd_img.o explore.o cudd_ckt.o
cudd_objs=$(CUDD)/obj/libobj.a $(CUDD)/cudd/libcudd.a $(CUDD)/mtr/libmtr.a $(CUDD)/st/libst.a $(CUDD)/util/libutil.a $(CUDD)/epd/libepd.a util/libcktutil.a

CUDDLIB_FLAGS=$(foreach libs,$(CUDDLIBS),-L$(CUDD)/$(libs)) $(foreach libs,$(CUDDLIBS),-l$(libs))

CUDDFLAGS:=-I$(CUDD)/include
CPPUTEST_FLAGS:=-I$(CPPUTEST_HOME)/include
CXXFLAGS:= $(CUDDFLAGS) $(CPPUTEST_FLAGS) -DCPU -g -Wall -std=c++11 -mtune=native -march=native -fopenmp -D_GLIBCXX_PARALLEL
LIBS:=-lstdc++ $(CUDDLIB_FLAGS) -L$(CPPUTEST_HOME) -L${CPPUTEST_HOME}/lib -L./util -lutil -lcudd -lobj -lmtr -lst -lepd -lcktutil -lm -lz -lCppUTestExt -lCppUTest 
CXX=g++
.PHONY: all

explore: $(objs) $(cudd_objs)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

explore.o: explore.cc cudd_ckt.h bdd_img.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
bdd_img.o: bdd_img.cc bdd_img.h cudd_ckt.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
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
	
test: AllTests.o cudd_ckt.o bdd_img.o getpis_test.o 
	@echo $(lflags)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $^
	./$@
