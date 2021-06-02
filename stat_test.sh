#!/bin/bash

./main -g oak_grove -stat_run 200 0.6  >>stat.txt 2>&1
./main -g oak_grove -stat_run 200 0.65 >>stat.txt 2>&1
./main -g oak_grove -stat_run 200 0.7  >>stat.txt 2>&1
./main -g oak_grove -stat_run 200 0.75 >>stat.txt 2>&1

./main -g oak_grove -stat_run 400 0.6  >>stat.txt 2>&1
./main -g oak_grove -stat_run 400 0.65 >>stat.txt 2>&1
./main -g oak_grove -stat_run 400 0.7  >>stat.txt 2>&1
./main -g oak_grove -stat_run 400 0.75 >>stat.txt 2>&1