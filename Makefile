CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp
LDFLAGS = -lhdf5_cpp -lhdf5

# Eigen library include path (adjust if different)
EIGEN_INC = /usr/include/eigen3

# Source and target
SRC = your_program.cpp
TARGET = your_program

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -I$(EIGEN_INC) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
