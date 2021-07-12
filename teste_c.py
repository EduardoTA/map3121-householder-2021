import AutovalsAutovecs as aa
import numpy as np
np.set_printoptions(linewidth=500)
np.set_printoptions(formatter={'float': '{: 12.10f}'.format})
import math as math

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)


f = 0
with open('input-c') as f:
    f = f.read().split()

n_total = int(f.pop(0)) # Número de nós da treliça
n_moveis = int(f.pop(0)) # Número de nós móveis da treliça
n_barras = int(f.pop(0)) # Número de barras da treliça

rho = float(f.pop(0)) # Densidade rho [kg/m³]
A = float(f.pop(0)) # Área transversal das barras [m²]
E = float(f.pop(0)) # Módulo de elasticidade em GPa
E = E*1e9 # Módulo de elasticidade em Pa

K = np.zeros((n_moveis*2, n_moveis*2)) # Matrix de rigidez total
m = np.zeros(n_moveis) # Vetor de massas dos nós
M = np.zeros((n_moveis*2, n_moveis*2))

count_nos = 0 # Contador para saber quando chegou nos nós fixos
for barra in range(0, n_barras):
    # Extraindo cada linha do arquivo correspondente às barras
    i = int(f[barra*4])
    i =i-1
    j = int(f[barra*4+1])
    j =j-1
    theta = float(f[barra*4+2])
    L = float(f[barra*4+3])

    alfa = A*E/L # Coeficiente da equação (1)
    C = math.cos(theta*math.pi/180)
    S = math.sin(theta*math.pi/180)

    print('i = '+str(i))
    print('j = '+str(j))
    print('theta = '+str(theta))
    print('L = '+str(L))
    print('C = '+str(C))
    print('S = '+str(S))
    


    # Fazendo o cálculo de cada elemento da matriz da equação (1)
    n = n_moveis*2
    if 0 <= 2*i < n:
        K[2*i, 2*i] = K[2*i, 2*i]+alfa*C**2
        if 0 <= 2*i+1 < n:
            K[2*i, 2*i+1] = K[2*i, 2*i+1]+alfa*C*S
            K[2*i+1, 2*i] = K[2*i+1, 2*i]+alfa*C*S
        if 0 <= 2*j < n:
            K[2*i, 2*j] = K[2*i, 2*j]-alfa*C**2
            K[2*j, 2*i] = K[2*j, 2*i]-alfa*C**2
        if 0 <= 2*j+1 < n:
            K[2*i, 2*j+1] = K[2*i, 2*j+1]-alfa*C*S
            K[2*j+1, 2*i] = K[2*j+1, 2*i]-alfa*C*S
    if 0 <= 2*i+1 < n:
        K[2*i+1, 2*i+1] = K[2*i+1, 2*i+1]+alfa*S**2
        if 0 <= 2*j < n:
            K[2*i+1, 2*j] = K[2*i+1, 2*j]-alfa*C*S
            K[2*j, 2*i+1] = K[2*j, 2*i+1]-alfa*C*S
        if 0 <= 2*j+1 < n:
            K[2*i+1, 2*j+1] = K[2*i+1, 2*j+1]-alfa*S**2
            K[2*j+1, 2*i+1] = K[2*j+1, 2*i+1]-alfa*S**2
    if 0 <= 2*j < n:
        K[2*j, 2*j] = K[2*j, 2*j]+alfa*C**2
        if 1 <= 2*j+1 < n:
            K[2*j, 2*j+1] = K[2*j, 2*j+1]+alfa*C*S
            K[2*j+1, 2*j] = K[2*j+1, 2*j]+alfa*C*S
    if 0 <= 2*j+1 < n:
        K[2*j+1, 2*j+1] = K[2*j+1, 2*j+1]+alfa*S**2
        
    try:
        m[i] = m[i]+ 1/2*rho*A*L
    except:
        pass
    try:
        m[j] = m[j]+ 1/2*rho*A*L
    except:
        pass
    
    print('\n')
    # Somando 1 no contador
    count_nos += count_nos

tilde_M = np.zeros((n_moveis*2, n_moveis*2))
for no in range(0, n_moveis):
    M[no*2,no*2] = m[no]
    M[no*2+1, no*2+1] = m[no]

    tilde_M[no*2,no*2] = 1/math.sqrt(m[no])
    tilde_M[no*2+1, no*2+1] = tilde_M[no*2,no*2]

tilde_K = np.matmul(np.matmul(tilde_M, K), tilde_M)

resultados = aa.AutovalsAutovecs(tilde_K, 1e-6, True)

for i in range(0, tilde_K.shape[0]):
    for j in range(0, tilde_K.shape[0]):
        if np.abs(tilde_K[i,j]) < 1e-6:
            tilde_K[i,j] = 0

print('massas')
print(m)

print(check_symmetric(M))

np.savetxt('K.txt', K)

np.savetxt('M.txt', M)

np.savetxt('tilde_M.txt', tilde_M)

# Calcula as frequências
w = np.zeros(n_moveis*2)
for i in range(n_moveis*2):
    w[i] = math.sqrt(resultados[0][i])

# Calcula os modos de vibração
z = np.zeros((n_moveis*2, n_moveis*2))
for j in range(n_moveis*2):
    z[:,j] = tilde_M.dot(resultados[1][:,j])

np.set_printoptions(precision=10)
for i in range(0, n_moveis*2):
    # Configuração de impressão para ter mais dígitos
    np.set_printoptions(formatter={'float': '{: 12.10f}'.format})

    # Impressão dos autovals/autovecs
    print('freq: {0:12.10f}, modo:'.format(w[i]), z[:, i])

print('M:')
print(M)
print('tilde_M:')
print(tilde_M)
print('tilde_K:')
print(tilde_K)

print(w)
