#!/bin/bash
make -j8 && ./main -sandbox
gprof ./main ./gmon.out > profile.txt
gprof2dot ./profile.txt | dot -Tsvg -o profile.svg