g++ -O3 -std=c++0x `root-config --cflags` interp.cpp `root-config --glibs` -lm -lnse -lgsl
chmod +x a.out
./a.out
