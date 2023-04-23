#!/bin/bash

cd ~/references_grade/differentiable-sdf-rendering/python

python3.8 optimize.py "$1" --optconfig diffuse-1 --force
python3.8 render_turntable.py "$1" --optconfig diffuse-1

python3.8 optimize.py "$1" --optconfig diffuse-2 --force
python3.8 render_turntable.py "$1" --optconfig diffuse-2

python3.8 optimize.py "$1" --optconfig diffuse-6 --force
python3.8 render_turntable.py "$1" --optconfig diffuse-6

python3.8 optimize.py "$1" --optconfig diffuse-12 --force
python3.8 render_turntable.py "$1" --optconfig diffuse-12

cd ~/grade