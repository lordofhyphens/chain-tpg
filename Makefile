CUDD:=./cudd
CUDDLIBS=cudd mtr util dddmp epd obj st 


CUDDFLAGS:=-I$(CUDD)/include   -lm
CXXFLAGS:= $(CUDDFLAGS) -DCPU -g -fopenmp -Wall -std=c++11 -mtune=native -L./util -lz
CXX=g++
.PHONY: all

explore: explore.cc cudd_ckt.h util/libcktutil.a
	$(CXX) $(CXXFLAGS) explore.cc -lcktutil $(CUDD)/obj/libobj.a $(CUDD)/cudd/libcudd.a $(CUDD)/mtr/libmtr.a $(CUDD)/st/libst.a $(CUDD)/util/libutil.a $(CUDD)/epd/libepd.a -o $@

clean: 
	rm explore
test: explore
	./$^ ../bench/s27.bench
	dot -Tpdf -O states.dot
	dot -Tpdf -O ckt.dot
	dot -Tpdf -O const.dot
