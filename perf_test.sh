#!/bin/bash
make -j8
#./main -demo 10 -no_debug -patch_size 10
make -j8 && ./main -sandbox -render
gprof ./main ./gmon.out > profile.txt
gprof2dot ./profile.txt | dot -Tsvg -o profile.svg