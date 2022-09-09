#NVFLAGS=--std=c++11 --gpu-architecture=sm_70 -O3
NVFLAGS=--std=c++11 -O3
ifeq ($(DEBUG_CPU),TRUE)
	NVFLAGS+=-DDEBUG_CPU
endif

ifeq ($(DEBUG_GPU),TRUE)
	NVFLAGS+=-DDEBUG_GPU
endif

ht_loc: main.o kernel.o helper.o
	hipcc $(NVFLAGS) -I. ./build/main.o ./build/helper.o ./build/kernel.o -o ./build/$@
main.o: main.cpp kernel.cpp helper.cpp
	hipcc $(NVFLAGS) -I. -c main.cpp -o ./build/$@ 
helper.o: helper.cpp helper.hpp
	hipcc $(NVFLAGS) -I. -c helper.cpp -o ./build/$@
kernel.o: kernel.cpp kernel.hpp
	hipcc $(NVFLAGS) -I. -c kernel.cpp -o ./build/$@
clean:
	rm ./build/*.o ./build/ht_loc
