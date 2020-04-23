#!/bin/bash
clear

gcc main3.c -o main

for file in  gr48
do
    for desv in 10 25 50
    do
        for inst in {1..10}
        do
            echo "./main < instances/$file"_"$desv"_"$inst.txt"
            ./main < instances/$file'_'$desv'_'$inst'.txt' >> Resultados2/"Resultados_$file"_"$desv"_"$inst.txt"
        done
    done
done