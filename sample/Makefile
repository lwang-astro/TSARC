## TARGET
TARGET=chain hierarchy

VPATH=../include

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -g -I../include -Wall
CXXFLAGS += -std=c++11

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## Time profile---------------------------------------#
CXXFLAGS += -D TIME_PROFILE

## for debug------------------------------------------#
#CXXFLAGS += -g -D DEBUG -fbounds-check -Wall

## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h extrapolation.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGET)
