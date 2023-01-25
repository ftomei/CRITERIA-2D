#!/bin/bash
cd src
python cythonSetup.py build_ext --inplace
cd ..
python src/main.py --path ./data/errano_evaluation_aug --iteration -1