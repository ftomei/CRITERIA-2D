#!/bin/bash
cd src
python cythonSetup.py build_ext --inplace
cd ..
python src/tuning.py --path ./data/errano_tuning --num_iterations 250 --seed 42