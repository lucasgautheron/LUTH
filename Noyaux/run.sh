g++ -O3 -std=c++0x `root-config --cflags` rates.cpp `root-config --glibs` -lm -lgsl -lgslcblas
chmod +x a.out
./a.out
