## TARGET
TARGET=chain hierarchy symtest keplerorbit

VPATH=../include

use_x86=no
use_k_computer=yes
QD_PATH=none
#QD_PATH=/opt/qd-2.3.18

## compiler-------------------------------------------#
ifeq ($(use_x86),yes)
CXX= g++
CXXFLAGS += -O2 -Wall
#CXXFLAGS += -g -fbounds-check
endif

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
CXXFLAGS += -Kfast
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
endif

## for debug------------------------------------------#
#CXXFLAGS += -D DEBUG 
#CXXFLAGS += -D ARC_DEBUG
#CXXFLAGS += -D ARC_DEEP_DEBUG

## Flag ----------------------------------------------#
CXXFLAGS += -I../include
CXXFLAGS += -std=c++11

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## profile---------------------------------------#
CXXFLAGS += -D ARC_PROFILE
## gperf
#GPERF_PATH=
#CXXLIBS += -L$(GPERF_PATH) -lprofiler -ltcmalloc

## qd-lib---------------------------------------------#
ifneq ($(QD_PATH),none)
CXXFLAGS += -D USE_QD ${shell ${QD_PATH}/bin/qd-config --cxxflags}
CXXLIBS += ${shell ${QD_PATH}/bin/qd-config --libs}
endif


## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h extrapolation.h Float.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

clean:
	rm -f $(TARGET)
