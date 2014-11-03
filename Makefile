CUDD:=./cudd
CUDDLIBS=cudd mtr util dddmp epd obj st 
objs=bdd_img.o explore.o cudd_ckt.o
cudd_objs=$(CUDD)/obj/libobj.a $(CUDD)/cudd/libcudd.a $(CUDD)/mtr/libmtr.a $(CUDD)/st/libst.a $(CUDD)/util/libutil.a $(CUDD)/epd/libepd.a util/libcktutil.a

CUDDFLAGS:=-I$(CUDD)/include
CXXFLAGS:= $(CUDDFLAGS) -DCPU -g -fopenmp -Wall -std=c++11 -mtune=native -march=native -D_GLIBCXX_PARALLEL
LIBS:=-L./util -lcktutil -lm -lz 
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
test: explore
	./$< ../bench/s27.bench
