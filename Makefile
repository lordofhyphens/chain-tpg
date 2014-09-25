CUDD:=./cudd
CUDDLIBS=cudd mtr util dddmp epd obj st 


CUDDFLAGS:=-I$(CUDD)/include   -lm
CXXFLAGS:= $(CUDDFLAGS) -g -fopenmp -Wall -std=c++11 -mtune=native -L./util
CXX=g++
.PHONY: all

explore:
	$(CXX) $(CXXFLAGS) explore.cc -lcktutil $(CUDD)/obj/libobj.a $(CUDD)/cudd/libcudd.a $(CUDD)/mtr/libmtr.a $(CUDD)/st/libst.a $(CUDD)/util/libutil.a $(CUDD)/epd/libepd.a 
