#!/bin/bash

for file in burma14 ulysses16 gr17 gr21 ulysses22 fri26 bayg29 gr48
do
    for desv in 10 25 50
    do
        for inst in {1..10}
        do
            echo "python < ResultadosAG_Artigo2/$file"_"$desv"_"$inst.txt"
            python3 ResultadosAG_Artigo2/"Resultados_$file"_"$desv"_"$inst"".csv" $file"_"$desv".csv"
        done
    done
done