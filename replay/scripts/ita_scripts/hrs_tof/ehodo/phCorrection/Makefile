# // Makefile //

OBJ=correction.o ParamMan.o
TARGET=phCorrection
INC=correction.hh ParamMan.hh

# for tfarm

#ROOT	= /cern/root/bin/root
#ROOTFLAGS = $(shell root-config --cflags)
#ROOTLIBS = $(shell root-config --libs)


ROOTFLAGS = -pthread -std=c++11 -m32 -msse -mfpmath=sse -I/home/dragon/sfw/root-6.08.06-install/./include
ROOTLIBS  = -L/home/dragon/sfw/root-6.08.06-install/./lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
#ROOTFLAGS = -pthread -std=c++11 -m32 -msse -mfpmath=sse -I/home/dragon/sfw/root-6.08.06-install/./include
#ROOTLIBS  = -L/home/dragon/sfw/root-6.08.06-install/./lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic                                                                     
CXX = /usr/local/bin/g++


#ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS) -lMinuit

#ROOT	= /home/dragon/sfw/root-6.08.06-install/bin/root
#ROOTFLAGS = $(shell ${ROOTSYS}/bin/root-config --cflags) 
#ROOTLIBS  = $(shell ${ROOTSYS}/bin/root-config --libs)
#CXXFLAGS = -Wall $(ROOTFLAGS)
#CXXLIBS  = $(ROOTLIBS) -lMinuit
#CXX = g++

all:$(TARGET)

$(TARGET):$(OBJ) $(INC)
	$(CXX) $(CXXLIBS) -o $(TARGET) $(OBJ)

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY:clean
clean:
	rm -f $(OBJ) $(TARGET) *~
