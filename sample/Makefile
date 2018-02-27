## TARGET
TARGET=chain hierarchy symtest keplerorbit

VPATH=../include

use_x86=yes
#use_k_computer=yes
QD_PATH=none
#QD_PATH=/opt/qd-2.3.18
GPERF_PATH=none
#GPERF_PATH = /opt/gperftools-2.6.1


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
CXXFLAGS += -DARC_OPT_SYM2

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## profile---------------------------------------#
CXXFLAGS += -D ARC_PROFILE

## gperf----------------------------------------------#
ifneq ($(GPERF_PATH),none)
CXXLIBS += -L$(GPERF_PATH)/lib -lprofiler -ltcmalloc
endif

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
