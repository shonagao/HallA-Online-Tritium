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

TARGET1= hrs_KID_F1_eff
OBJS1=   hrs_KID_F1_eff.o  
TARGET2= hrs_KID_F1
OBJS2=   hrs_KID_F1.o  
TARGET3= hrs_path_corr
OBJS3=   hrs_path_corr.o  
TARGET4= hrs_KID_F1_eff2
OBJS4=   hrs_KID_F1_eff2.o  
TARGET5= hrs_KID_F1_check
OBJS5=   hrs_KID_F1_check.o  
TARGET6= hrs_MM_F1
OBJS6=   hrs_MM_F1.o  

.PHONY: all
all: $(TARGET1)\
     $(TARGET2)\
     $(TARGET3)\
     $(TARGET4)\
     $(TARGET5)\
     $(TARGET6)
$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*.d core $(BINDIR)/*

