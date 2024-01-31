#!/bin/bash

cd ~/references_grade/differentiable-sdf-rendering/python

cp ~/grade_resources/Huawei_models/h1_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_1

cp ~/grade_resources/Huawei_models/h2_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_2

cp ~/grade_resources/Huawei_models/h3_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_3

cp ~/grade_resources/Huawei_models/h4_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_4

cp ~/grade_resources/Huawei_models/h5_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_5

cp ~/grade_resources/Huawei_models/h6_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_6

cp ~/grade_resources/Huawei_models/h7_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_7

cp ~/grade_resources/Huawei_models/h8_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_8

cp ~/grade_resources/Huawei_models/h9_norm.obj ~/references_grade/differentiable-sdf-rendering/scenes/d1/d1.obj
python optimize.py "d1" --optconfig diffuse-12 --force
python render_turntable.py "d1" --optconfig diffuse-12
cp -r ~/references_grade/differentiable-sdf-rendering/outputs/d1 ~/references_grade/differentiable-sdf-rendering/outputs/detail_9

cd ~/grade