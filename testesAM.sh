#!/bin/bash
clear

gcc mainAM.c -o mainAM

for file in burma14 ulysses16 gr17 gr21 ulysses22 fri26 bayg29 gr48
do
    for desv in 10 25 50
    do
        for inst in {1..10}
        do
            echo "./mainAM < instances/$file"_"$desv"_"$inst.txt"
            ./mainAM < instances/$file'_'$desv'_'$inst'.txt' > ResultadosAM_Artigo2/"Resultados_$file"_"$desv"_"$inst.csv"
        done
    done
done