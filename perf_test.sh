#!/bin/bash
make -j8
./main -demo 100 -no_debug -patch_size 10
gprof ./main ./gmon.out >> profile.txt
gprof2dot ./profile.txt | dot -Tsvg -o profile.svg