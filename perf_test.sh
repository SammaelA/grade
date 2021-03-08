#!/bin/bash

make && ./main -g grove160 -only_gen >>t160.txt 2>&1 && gprof ./main ./gmon.out >> t160.txt
make && ./main -g grove320 -only_gen >>t320.txt 2>&1 && gprof ./main ./gmon.out >> t320.txt