#!/bin/bash
make && ./main -g grove10 -only_gen >s10.txt 2>&1
make && ./main -g grove50 -only_gen >s50.txt 2>&1
make && ./main -g grove100 -only_gen >s100.txt 2>&1
make && ./main -g grove150 -only_gen >s150.txt 2>&1
make && ./main -g grove200 -only_gen >s200.txt 2>&1
make && ./main -g grove250 -only_gen >s250.txt 2>&1