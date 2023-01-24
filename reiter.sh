#!/bin/bash
until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_6.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done
