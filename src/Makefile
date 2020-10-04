ht_loc: main.o kernel.o helper.o
	nvcc --gpu-architecture=sm_61 -I. ./build/main.o ./build/helper.o ./build/kernel.o -o ./build/$@
main.o: main.cpp kernel.cpp helper.cpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c main.cpp -o ./build/$@ 
helper.o: helper.cpp helper.hpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c helper.cpp -o ./build/$@
kernel.o: kernel.cpp kernel.hpp
	nvcc -x cu --gpu-architecture=sm_61 -I. -c kernel.cpp -o ./build/$@
clean:
	rm ./build/*.o ./build/ht_loc