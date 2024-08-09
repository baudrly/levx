CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -march=native -I/usr/include/eigen3/
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5

SRC = your_program.cpp
TARGET = your_program

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -I$(EIGEN_INC) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
