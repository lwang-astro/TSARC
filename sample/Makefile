## TARGET
TARGET=chain hierarchy

VPATH=../include

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -I../include -Wall
CXXFLAGS += -std=c++11

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## Time profile---------------------------------------#
CXXFLAGS += -D TIME_PROFILE

## for debug------------------------------------------#
CXXFLAGS += -O2 
#CXXFLAGS += -g -fbounds-check
#CXXFLAGS += -D DEBUG 
CXXFLAGS += -D ARC_DEBUG

## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h extrapolation.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGET)
