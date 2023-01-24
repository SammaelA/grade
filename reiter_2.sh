#!/bin/bash
until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_1.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_2.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_3.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_5.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/cup_6.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/plate_9.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done

until ./main -sandbox -opt_benchmark diff_optimization_benchmark.blk resources/textures/DishesData/plate_10.jpg
do
    if test -f "config/backup.blk"
    then
        echo "FALL"
    else
        break
    fi
done