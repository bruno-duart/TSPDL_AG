import pandas as pd 
import sys

dataF = pd.read_csv(f'{sys.argv[1]}.csv')
#dataF = pd.read_csv(f'ResultadosAG_Artigo2/Resultados_burma14_50_2.csv')
#print(dataF)
df = dataF.nlargest(11,'MediaRes').tail(1)
print(df)



#.query(f'MediaRes == {dataF["MediaRes"].min()}')