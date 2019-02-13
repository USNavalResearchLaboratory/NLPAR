
#CSOURCES =  djr_npar.c #cNLPAR.c 
CPPSOURCES= cNLPARv5.cpp 

OBJECTS= $(CSOURCES:.c=.o) $(CPPSOURCES:.cpp=.o)
EXECUTABLE = patternProcessingC.dylib

CFLAGS= -O3 -m64 -c -ffast-math -ftree-vectorize -march=native -static-libgcc 
CPPFLAGS= -O3 -m64 -c -std=c++11 -ffast-math -ftree-vectorize -march=native -static-libstdc++


LFLAGS=  -m64 -O3 --shared 
#LFLAGS=  -m64 --shared 
#FTNCOMP = gfortran

OPENMP=-Xpreprocessor -fopenmp -lomp
CCOMP= clang 
CPPCOMP= clang++
#OPENMP = -fopenmp
#CCOMP = /usr/local/bin/gcc
#CPPCOMP = /usr/local/bin/g++

all: $(CSOURCES) $(CPPSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPPCOMP) $(LFLAGS) $(OPENMP) $(OBJECTS) -o $@

.c.o:
	$(CCOMP) $(CFLAGS) $(OPENMP) $< -o $@

.cpp.o:
	$(CPPCOMP) $(CPPFLAGS) $(OPENMP) $< -o $@
	
clean:
	rm -f $(OBJECTS)

# project(patternProcessing)
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp" CACHE STRING "" FORCE)
# set(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp" CACHE STRING "" FORCE)
# add_library(patternProcessing SHARED djr_npar.c cNLPAR.c Cpp_Cexample.cpp)
