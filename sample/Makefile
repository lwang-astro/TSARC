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
## gperf
#GPERF_PATH=
#CXXLIBS += -L$(GPERF_PATH) -lprofiler -ltcmalloc

## qd-lib---------------------------------------------#
QD_PATH=/opt/qd-2.3.18
CXXFLAGS += -D USE_QD ${shell ${QD_PATH}/bin/qd-config --cxxflags}
CXXLIBS += ${shell ${QD_PATH}/bin/qd-config --libs}

## for debug------------------------------------------#
CXXFLAGS += -O2 
#CXXFLAGS += -g -fbounds-check
#CXXFLAGS += -D DEBUG 
#CXXFLAGS += -D ARC_DEBUG
#CXXFLAGS += -D ARC_DEEP_DEBUG

## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h extrapolation.h Float.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

clean:
	rm -f $(TARGET)
