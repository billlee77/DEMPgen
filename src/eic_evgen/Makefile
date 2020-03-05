# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
# 
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --cflags --glibs)
ROOTINC       = $(shell root-config --incdir)
#CXXFLAGS      = -Wall -g -frtti -fexceptions -fPIC -O 
CXXFLAGS      = -Wall -frtti -fexceptions -fPIC -O 
O=eic

# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include -I$(ROOTINC) 
CXX           = clang++
#CXX           = g++
LD            = clang++
#LD            = g++
LDFLAGS       = 

LIBS          = $(ROOTLIBS) 
GLIBS         = $(ROOTGLIBS)

ALL_LIBS =  $(GLIBS) $(LIBS)

# Test code executibles
PROGS = $(O)

$(O): $(O).cxx
	rm -f $@
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@  $(O).cxx $(ALL_LIBS)

         
clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###
