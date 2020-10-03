ht_loc: main.cpp kernel.* helper.*
	nvcc -x cu --gpu-architecture=sm_61 main.cpp helper.cpp kernel.cpp -o $@