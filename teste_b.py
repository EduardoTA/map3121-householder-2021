import AutovalsAutovecs as aa
import numpy as np
import math as math
from teste_a import EstimativasErroGerais as EstimativasErroGerais

def EstimativasErroAnalitico(Lambda):
    print("===Estimativas de Erro Analíticas===\n")

    n = Lambda.shape[0]

    # Calcula os autovalores reais, pela fórmula analítica, e ordena eles do maior para o menor
    autovalores_reais = 1/2*(1-np.cos((2*np.arange(n, 0, -1)-1)*np.pi/(2*n+1)))**(-1)
    autovalores_reais = np.flip(np.sort(autovalores_reais))

    autovalores_obtidos = np.copy(np.flip(np.sort(Lambda)))

    erros = np.zeros(n)

    # Faz a comparação entre os autovalores reais e os obtidos pelo método
    print('Autovalor obtido | Autovalor real | erro')
    for i in range(0, n):
        erros[i] = math.sqrt(math.pow(autovalores_obtidos[i] - autovalores_reais[i], 2))
        print('{0:12.10f}       {1:12.10f}     {2}'.format(autovalores_obtidos[i], autovalores_reais[i],erros[i]))
    print('-----------------------------')
    print('Erro máx: ', np.max(erros))
    print('\n')


def teste_b():
    # Menu de seleção
    print('\n############')
    print('Teste b selecionado')
    print('############\n')

    epsilon = 0.1
    deslocamentos = True
    n = 20

    try:
        epsilon = float(input('epsilon = '))
        n = int(input('n = '))
        deslocamentos = input('Usar deslocamentos espectrais? (s,n)')
        if deslocamentos == 's':
            deslocamentos = True
        else:
            deslocamentos = False
    except:
        print('epsilon deve ser número real, pe: 1e-6 e n deve ser inteiro')
        return

    print("\n======================RESULTADOS DO TESTE B======================")

    A = np.ones((n,n))
    for k in range(0,n):
        for j in range(0,k+1):
            A[k,j] = n-k
        for i in range(0,k):
            A[i,k] = n-k

    resultados = aa.AutovalsAutovecs(A, epsilon, deslocamentos)

    print("\nForam necessárias k = {0} iterações no método QR\n".format(resultados[2]))

    for i in range(0, A.shape[0]):
        # Configuração de impressão para ter mais dígitos
        #np.set_printoptions(formatter={'float': '{: 12.10ff}'.format})

        # Impressão dos autovals/autovecs
        print('Autovalor: {0:12.10f}, Autovetor:'.format(resultados[0][i]), resultados[1][:, i])
    
    print("\n")
    EstimativasErroGerais(A, resultados[0], resultados[1])
    EstimativasErroAnalitico(resultados[0])
