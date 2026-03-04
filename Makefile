CXX      := g++
CXXFLAGS := -O2 -g -Wall -Wextra -std=c++17 -fPIC

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

INC := -Iinclude \
       -I$(WCSIM_BUILD_DIR)/include \
       -I$(WCSIM_BUILD_DIR)/include/WCSim \
       -I$(WCSIM_BUILD_DIR)/include/WCSimRoot \
       -I$(WCSIM_SOURCE_DIR)/include \
       -I$(WCSIM_SOURCE_DIR)/include/WCSim \
       -I$(WCSIM_SOURCE_DIR)/include/WCSimRoot

WCSIM_LIB := -L$(WCSIM_BUILD_DIR)/lib -lWCSimRoot
RPATH     := -Wl,-rpath,$(WCSIM_BUILD_DIR)/lib

TARGET := beam_mc_studies

SRCS := main.cc \
        src/Run.cc \
        src/Utils.cc \
        src/Studies_BasicSpectra.cc \
        src/Studies_QperHit.cc \
        src/Studies_QvsN_2D.cc \
        src/Studies_VtxXZ.cc \
		src/Studies_Cone42.cc

OBJS := $(SRCS:.cc=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(ROOTLIBS) $(WCSIM_LIB) $(RPATH) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INC) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)