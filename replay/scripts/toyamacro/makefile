CC = gcc 
CXX = g++ 
CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1=     s2_t0_calib
OBJS1=       s2_t0_calib.o Tree.o Setting.o ParamMan.o
TARGET2=     RHRS_PID
OBJS2=       RHRS_PID.o Tree.o Setting.o ParamMan.o
TARGET3=     s2_twc
OBJS3=       s2_twc.o Tree.o Setting.o ParamMan.o
TARGET4=     mk_twc_tree
OBJS4=       mk_twc_tree.o Tree.o ParamMan.o
TARGET5=     draw_s2tof
OBJS5=       draw_s2tof.o Tree.o ParamMan.o Setting.o
TARGET6=     fbus_t0_rough
OBJS6=       fbus_t0_rough.o Tree.o ParamMan.o Setting.o
TARGET7=     ana_cointime
OBJS7=       ana_cointime.o Tree.o Setting.o ParamMan.o
TARGET8=     calc_charge
OBJS8=       calc_charge.o
TARGET9=     f1_t0_rough
OBJS9=       f1_t0_rough.o Tree.o ParamMan.o Setting.o
TARGET10=    ana_cointime_ac1study
OBJS10=      ana_cointime_ac1study.o Tree.o ParamMan.o Setting.o

all: $(TARGET1) \
     $(TARGET2) \
     $(TARGET3) \
     $(TARGET4) \
     $(TARGET5) \
     $(TARGET6) \
     $(TARGET7) \
     $(TARGET8) \
     $(TARGET9) \
     $(TARGET10) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

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

$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

$(TARGET8): $(patsubst %,$(LIBDIR)/%,$(OBJS8))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

$(TARGET9): $(patsubst %,$(LIBDIR)/%,$(OBJS9))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

$(TARGET10): $(patsubst %,$(LIBDIR)/%,$(OBJS10))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

