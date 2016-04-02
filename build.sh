#!/bin/bash

#source ./version.sh
git pull
make clean
time make -j2
./src/xmolview
