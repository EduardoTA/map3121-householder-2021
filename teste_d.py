import AutovalsAutovecs as aa
import numpy as np
import math as math
from teste_a import EstimativasErroGerais as EstimativasErroGerais

def teste_d():
    # Menu de seleção
    print('\n############')
    print('Teste d selecionado')
    print('############\n')

    epsilon = 0.1
    deslocamentos = True
    f = 0
    n = 0
    A = np.zeros((n,n))

    try:
        nome_arquivo = str(input('Digite o nome do arquivo: '))
        with open(nome_arquivo) as f:
            f = f.read().split()
    except:
        print('Arquivo não existe ou está indisponível')
        return

    try:
        epsilon = float(input('epsilon = '))
        deslocamentos = input('Usar deslocamentos espectrais? (s,n)')
        if deslocamentos == 's':
            deslocamentos = True
        else:
            deslocamentos = False
    except:
        print('epsilon deve ser número real, pe: 1e-6')
        return

    try:
        n = int(f.pop(0)) # Número de linhas da matriz quadrada
        A = np.zeros((n,n))
    except:
        print('Erro na leitura do número de linhas da matriz')
        return
    try:
        for i in range(0,n):
            for j in range(0,n):
                A[i,j] = float(f[j+i*n])
    except:
        print('Arquivo mal formatado')
        return
    
    print("\n======================RESULTADOS DO TESTE D======================")

    resultados = aa.AutovalsAutovecs(A, epsilon, deslocamentos)

    print("\nForam necessárias k = {0} iterações no método QR\n".format(resultados[2]))

    for i in range(0, A.shape[0]):
        # Configuração de impressão para ter mais dígitos
        #np.set_printoptions(formatter={'float': '{: 12.10f}'.format})

        # Impressão dos autovals/autovecs
        print('Autovalor: {0:12.10f}, Autovetor:'.format(resultados[0][i]), resultados[1][:, i])
    
    print("\n")
    EstimativasErroGerais(A, resultados[0], resultados[1])

    #Lambda = np.zeros((A.shape[0], A.shape[0]))
    #for i in range(0, A.shape[0]):
    #    Lambda[i,i] = resultados[0][i]
    
    #erro = np.max(np.abs((np.matmul(np.matmul(resultados[1], Lambda), np.transpose(resultados[1])) - A)))
    #print("\nEstimativa de erro H*L*HT-A = {0}\n".format(erro))

    return



