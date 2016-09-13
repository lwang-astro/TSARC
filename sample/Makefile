## TARGET
TARGET=chain

## compiler-------------------------------------------#
CXX= g++

## Flag ----------------------------------------------#
CXXFLAGS = -O2 -I../include

## for debug------------------------------------------#
CXXFLAGS += -g -D DEBUG -fbounds-check -Wall

## Target --------------------------------------------#

$(TARGET): %: %.cxx
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGET)
