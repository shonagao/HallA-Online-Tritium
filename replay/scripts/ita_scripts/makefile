CC=gcc
CXX=g++

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

# mp (Make Parameter)
#TARGET1= hrs_tof_twcorr
#OBJS1=   hrs_tof_twcorr.o 

TARGET1= conv
OBJS1=   conv.o  
TARGET2= conv_TOF
OBJS2=   conv_TOF.o  
TARGET3= copy
OBJS3=   copy.o  

.PHONY: all
all: $(TARGET1)\
     $(TARGET2)\
     $(TARGET3)
$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*.d core $(BINDIR)/*

