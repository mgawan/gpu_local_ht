ht_loc: main.cpp
	nvcc -x cu main.cpp helper.cpp kernel.cpp -o $@