## TARGET
TARGET=chain hierarchy symtest

VPATH=../include

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -I../include -Wall
CXXFLAGS += -std=c++11

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## profile---------------------------------------#
CXXFLAGS += -D ARC_PROFILE
#CXXLIBS += -L$(GPERF_PATH) -lprofiler -ltcmalloc

## for debug------------------------------------------#
CXXFLAGS += -O2 
#CXXFLAGS += -g -fbounds-check
#CXXFLAGS += -D DEBUG 
#CXXFLAGS += -D ARC_DEBUG

## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h extrapolation.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

clean:
	rm -f $(TARGET)
