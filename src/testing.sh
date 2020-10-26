#!/bin/bash

srun ./build/ht_loc ../locassm_data/localassm_extend_7-21.dat

diff_test=$(diff ./contig-test.dat ../locassm_data/res_localassm_extend_7-21.dat)

if ["$diff_test" == ""]
then
    echo "TEST PASSED!!"
    exit 0
else
    echo "TEST FAILED!!, check log for results"
    echo "$diff_test" >> log_test
    exit 1
fi