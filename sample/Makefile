## TARGET
TARGET=chain hierarchy

VPATH=../include

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -O2 -I../include

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
