
FLAGS=-fPIC -std=c++11

default: all

dirs:
	mkdir -p bin

spherical_harmonics: dirs sh/spherical_harmonics.cc
	$(CXX) $(FLAGS) -c sh/spherical_harmonics.cc -I./sh -o bin/spherical_harmonics.o
	$(CXX) $(FLAGS) -shared sh/spherical_harmonics.cc -I./sh -o bin/spherical_harmonics.so
	ar rcs bin/spherical_harmonics.a bin/spherical_harmonics.o

all: spherical_harmonics




