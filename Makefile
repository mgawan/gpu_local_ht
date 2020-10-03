ht_loc: main.o kernel.o helper.o
	nvcc --gpu-architecture=sm_61 -I. main.o helper.o kernel.o -o $@
main.o: main.cpp kernel.cpp helper.cpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c main.cpp -o $@ 
helper.o: helper.cpp helper.hpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c helper.cpp -o $@
kernel.o: kernel.cpp kernel.hpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c kernel.cpp -o $@
clean:
	rm *.o