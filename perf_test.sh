#!/bin/bash
make && ./main -g grove10 -only_gen >>t10.txt 2>&1 && gprof ./main ./gmon.out >> t10.txt
make && ./main -g grove20 -only_gen >>t20.txt 2>&1 && gprof ./main ./gmon.out >> t20.txt
make && ./main -g grove40 -only_gen >>t40.txt 2>&1 && gprof ./main ./gmon.out >> t40.txt
make && ./main -g grove80 -only_gen >>t80.txt 2>&1 && gprof ./main ./gmon.out >> t80.txt
make && ./main -g grove160 -only_gen >>t160.txt 2>&1 && gprof ./main ./gmon.out >> t160.txt
make && ./main -g grove320 -only_gen >>t320.txt 2>&1 && gprof ./main ./gmon.out >> t320.txt