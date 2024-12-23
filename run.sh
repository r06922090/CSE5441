#!/bin/bash

LANG=en_US.utf8
LC_ALL=en_US.utf8

REAL="real.exe"
mpicc -O2 cannon.cpp -o $REAL -lm

mpirun -np 1 ./$REAL
mpirun -np 16 ./$REAL
