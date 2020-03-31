#!/bin/bash

INSTANCE=gr21

echo -e "  "
gcc main2.c -o main -DMAX_ITER=10
./main < instances/"$INSTANCE"_50_1.txt

echo -e "  "
gcc main2.c -o main -DMAX_ITER=1000
./main < instances/"$INSTANCE"_50_1.txt