
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

dados = {"TaxaMut":[] , "MediaRes":[], "PercMenor":[], "MediaCont": []}
with open("Resultados_Burma14_10_1.txt", "r") as f:
    for line in f:
        line = line.split()
        dados["TaxaMut"].append(line[1])
        dados["MediaRes"].append(line[6])
        dados["MediaCont"].append(line[9])
        dados["PercMenor"].append(line[14])

TaxaMut = np.array(dados["TaxaMut"])
MediaRes = np.array(dados["MediaRes"])
PercMenor = np.array(dados["PercMenor"])
MediaCont = np.array(dados["MediaCont"])

plt.figure
plt.bar(TaxaMut, MediaRes)#MediaRes)
plt.title("Média dos resultados")
plt.axis([1, 20, min(MediaRes), max(MediaRes)])
plt.show()
plt.bar(TaxaMut, PercMenor)
plt.title("Percentual de vezes que o menor resultado aparece")
plt.show()
plt.bar(TaxaMut, MediaCont)
plt.title("Número médio de iterações")
plt.show()

