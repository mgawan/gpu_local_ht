NVFLAGS=--std=c++11

ht_loc: main.o kernel.o helper.o
	nvcc --gpu-architecture=sm_70 $(NVFLAGS) -I. ./build/main.o ./build/helper.o ./build/kernel.o -o ./build/$@
main.o: main.cpp kernel.cpp helper.cpp
	nvcc -x cu --gpu-architecture=sm_70 $(NVFLAGS) -I. -c main.cpp -o ./build/$@ 
helper.o: helper.cpp helper.hpp
	nvcc -x cu --gpu-architecture=sm_70 $(NVFLAGS) -I. -c helper.cpp -o ./build/$@
kernel.o: kernel.cpp kernel.hpp
	nvcc -x cu --gpu-architecture=sm_70 $(NVFLAGS) -I. -c kernel.cpp -o ./build/$@
clean:
	rm ./build/*.o ./build/ht_loc