#!/bin/bash
clear

gcc mainAG.c -o mainAG

for file in burma14 ulysses16 gr17 gr21 ulysses22 fri26 bayg29 gr48
do
    for desv in 10 25 50
    do
        for inst in {1..10}
        do
            echo "./mainAG < instances/$file"_"$desv"_"$inst.txt"
            time ./mainAG < instances/$file'_'$desv'_'$inst'.txt' > ResultadosAG_Artigo2_Ordenado/"Resultados_$file"_"$desv"_"$inst.csv"
        done
    done
done