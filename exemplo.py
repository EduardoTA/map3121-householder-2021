# Importa o package com a função que calcula os autovalores e autovetores
import AutovalsAutovecs as aa

import numpy as np
import math as math

# Matriz de exemplo
A = np.matrix([[2,4,1,1],[4,2,1,1],[1,1,1,2],[1,1,2,1.0]])

# AutovalsAutovecs é o método que calcula os autovalores e autovetores da matriz simétrica A
# Args.:
# A: Matriz simétrica (np.matrix)
# epsilon: tolerância (float)
# deslocamentos: usar deslocamentos espectrais (bool)
#
# Retorna uma tupla cujas posições são:
# resultados[0]: Vetor com os autovalores (np.array)
# resultados[1]: Matriz cujas colunas são os autovetores de A (np.matrix)
resultados = aa.AutovalsAutovecs(A, 1e-6, True)
