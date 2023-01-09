#!/bin/bash
until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk
do
echo "FALL\n"
done
fi
