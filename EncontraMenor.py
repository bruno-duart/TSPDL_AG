import pandas as pd 
import sys

instances = ['bayg29','burma14','fri26','gr17','gr21','ulysses16','ulysses22','gr48']
std = [10, 25, 50]

for inst in instances:
    for std_i in std:
        for i in range(1,11):
            dataF = pd.read_csv(f'ResultadosAG_Artigo2/Resultados_{inst}_{std_i}_{i}.csv')
            #dataF = pd.read_csv(f'ResultadosAG_Artigo2/Resultados_burma14_50_2.csv')
            #print(dataF)
            df = dataF.nlargest(11,'MediaRes').tail(1)
           # print(df)
            with open(f'ResultadosAG_Artigo2/ResultadosMédia_{inst}_{std_i}.csv', 'a') as f:
                df.to_csv(f, header=False)



#.query(f'MediaRes == {dataF["MediaRes"].min()}')