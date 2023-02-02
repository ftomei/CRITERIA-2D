#!/bin/bash
cd src
python cythonSetup.py build_ext --inplace
cd ..
python src/wrapper.py --path ./data/$1 --iteration $2 --seed 42