## TARGET
TARGET=chain

VPATH=../include

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -O2 -I../include

## for debug------------------------------------------#
#CXXFLAGS += -g -D DEBUG -fbounds-check -Wall

## Target --------------------------------------------#

$(TARGET): %: %.cxx AR.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGET)
