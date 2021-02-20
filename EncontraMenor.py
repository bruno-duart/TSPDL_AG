import pandas as pd 
import sys

dataF = pd.read_csv(f'{sys.argv[1]}.csv')
print(dataF)

print(dataF['MediaRes'].min())
