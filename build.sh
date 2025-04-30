#!/bin/bash

DEFAULTTHREADS=18
THREADS=${1:-$DEFAULTTHREADS}

echo "Using $THREADS threads"
cd /mnt/blocksci

(mkdir -p build && \
    cd build && \
    CC=gcc-7 CXX=g++-7 cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$THREADS && \
    make install) || exit 1

cd /mnt/blocksci

pip3 install -r pip-all-requirements.txt || exit 1

pip3 install jupyter notebook
pip3 install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user


CC=gcc-7 CXX=g++-7 pip3 install -e blockscipy || exit 1

cd Notebooks

jupyter notebook --ip="0.0.0.0" --allow-root || exit 1