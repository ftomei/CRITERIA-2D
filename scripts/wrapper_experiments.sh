#!/bin/bash
cd src
python cythonSetup.py build_ext --inplace
cd ..
python src/wrapper.py -p data/$1 -it $2