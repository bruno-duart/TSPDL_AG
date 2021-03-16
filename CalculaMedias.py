# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 21:41:40 2021

@author: bruno
"""

import pandas as pd

instances = ['bayg29','burma14','fri26','gr17','gr21',
             'ulysses16','ulysses22','gr48']
std = [10, 25, 50]

MediaRes = []
T_Med = []
Err_Med = []
BestRes = []
T_Best = []
Err_Best = []

for inst in instances:
    for std_i in std:
        with open(f'ResultadosAM_Artigo2/ResultadosMédia_{inst}_{std_i}.csv', 'r') as f:
            df = pd.read_csv(f)
            MediaRes.append(df['MediaRes'].mean())
            T_Med.append(df['TempoMed'].mean())
            Err_Med.append(df['ErroMed'].mean())
            BestRes.append(df['BestRes'].mean())
            T_Best.append(df['TempoMenor'].mean())
            Err_Best.append(df['ErroMenor'].mean())

df_d = {'Instância':['bayg29_10', 'bayg29_25', 'bayg29_50',
                     'burma14_10', 'burma14_25', 'burma14_50',
                     'fri26_10', 'fri26_25', 'fri26_50',
                     'gr17_10', 'gr17_25', 'gr17_50',
                     'gr21_10', 'gr21_25', 'gr21_50',
                     'ulysses16_10', 'ulysses16_25', 'ulysses16_50',
                     'ulysses22_10', 'ulysses22_25', 'ulysses22_50',
                     'gr48_10', 'gr48_25', 'gr48_50'], 
        'Sol. Méd.': MediaRes, 'T_med(s)': T_Med, '$\Delta_{Med} (\%)$': Err_Med,
        'Sol. Mel.': BestRes, 'T_mel(s)': T_Best, '$\Delta_{Me} (\%)$': Err_Best}
DF_result = pd.DataFrame(df_d)

DF_result.to_csv('ResultadoAM_Art2.csv',index=False)
print('Acabou')