#!/bin/bash
clear

gcc main3.c -o main

for file in burma14 ulysses16 gr17 gr21 ulysses22 fri26 bayg29 gr48
do
    for desv in 10 25 50
    do
        for inst in {1..10}
        do
            echo "./main < instances/$file"_"$desv"_"$inst.txt"
            ./main < instances/$file'_'$desv'_'$inst'.txt' >> Resultados/"Resultados_$file"_"$desv"_"$inst.txt"
        done
    done
done