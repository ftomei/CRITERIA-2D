#!/bin/bash
cd src
python cythonSetup.py build_ext --inplace
cd ..
python src/wrapper.py --path ./data/errano_tuning --num_iterations 500 --seed 42