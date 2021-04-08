#!/bin/bash
for (( i=1; i <= 1000; i++ ))
do
echo $i
nvidia-smi
sleep 1
done
