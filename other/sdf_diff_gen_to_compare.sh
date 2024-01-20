#!/bin/bash

cd ~/references_grade/differentiable-sdf-rendering/python

python optimize.py "$1" --optconfig diffuse-1 --force
python render_turntable.py "$1" --optconfig diffuse-1

python optimize.py "$1" --optconfig diffuse-2 --force
python render_turntable.py "$1" --optconfig diffuse-2

python optimize.py "$1" --optconfig diffuse-6 --force
python render_turntable.py "$1" --optconfig diffuse-6

python optimize.py "$1" --optconfig diffuse-12 --force
python render_turntable.py "$1" --optconfig diffuse-12

cd ~/grade