ht_loc: main.cpp kernel.* helper.*
	nvcc -x cu main.cpp helper.cpp kernel.cpp -o $@