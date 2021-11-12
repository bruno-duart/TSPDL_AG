"""
SCRIPT DE MIGRAÇÃO TXT -> CSV

Este script migra os resultados armazenados em txt, dos experimentos antigos,
para o formato csv, compatível com os novos experimentos.

Assim, também, facilita o tratamento dos dados e análises estatísticas.

"""

# libraries
import os
import numpy as np
import pandas as pd

# presets
SOURCE_PATH = 'Resultados2/'
TARGET_PATH = 'ResultadosAM_Artigo2/'
HEADER = {
    'Mutacao': [],
    'MediaRes': [],
    'NumMedGer': [],
    'BestRes': [],
    'ErroMed': [],
    'ErroMenor': [],
    'TempoMed': [],
    'TempoMenor': []
}

# 1) ler todos os arquivos de uma pasta
os.chdir('../'+SOURCE_PATH)
for filename in os.listdir():
	name = filename.split('.')[0]
	print(name)
	for k in HEADER.keys():
	    HEADER[k] = []
    # 2) para cada arquivo, iterar as linhas
	with open(filename) as f:
		lines = f.readlines()[3:]
		for ln in lines:
			mutacao = int(ln[1:5].replace(' ',''))
			ln = ln.replace('(',' ')
			values = [ s for s in ln[5:].split(' ') if s not in [' ',''] ]
			HEADER['Mutacao'].append( mutacao )
			HEADER['MediaRes'].append( int(values[4]) )
			HEADER['NumMedGer'].append( int(values[7]) )
			HEADER['BestRes'].append( int(values[10]) )
			HEADER['ErroMed'].append( int(values[11]) )
			HEADER['ErroMenor'].append( float(values[15]) )
			HEADER['TempoMed'].append( float(values[18]) )
			HEADER['TempoMenor'].append( float(values[21]) )
	# 3) exportar um csv
	df = pd.DataFrame(HEADER)
	df.to_csv('../'+TARGET_PATH+name+'.csv', index=False)
