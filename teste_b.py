import AutovalsAutovecs as aa
import numpy as np
import math as math

def EstimativasErroGerais(A, Lambda, HV):
    n = Lambda.shape[0]

    # Faz uma estimativa de erro = max(abs(H*V*VT*HT-I))
    erro = np.max(np.abs(np.matmul(HV, np.transpose(HV)) - np.identity(n)))
    print('max(abs(H*V*VT*HT-I)) = {0}'.format(erro))

    # Faz uma estimativa de erro = max(abs(matmul(A,H*V)-matmul(H*V,L)))
    L = np.zeros((n, n))

    for i in range(0, n):
        L[i, i] = Lambda[i]

    erro = np.max(np.abs(np.matmul(A, HV) - np.matmul(HV, L)))
    print('max(abs(matmul(A,H*V)-matmul(H*V,L))) = {0}\n'.format(erro))

    # Faz uma estimativa de erro comparando (A*v)/v, com os autovalores obtidos
    erros = np.zeros(n)
    for i in range(0, n):
        erros[i] = math.sqrt(math.pow(np.mean(np.divide(np.matmul(A, HV[:, i]), HV[:, i])) - Lambda[i], 2))
        print('Autovalor obtido fazendo matdiv(matmul(A,v),v): {0:12.10f},'
                '  Autovalor obtido pelo método: {1:12.10f}, erro: {2}'
                .format(np.mean(np.divide(np.matmul(A, HV[:, i]), HV[:, i])), Lambda[i], erros[i]))
    print('Erro máximo = {0}\n'.format(np.max(erros)))

def EstimativasErroAnalitico(Lambda):
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
    print('############')
    print('Teste b selecionado')
    print('############')

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

    print("\nOs autovalores e autovetores da matriz do item 4.1.b serão impressos\n")

    A = np.ones((n,n))
    for k in range(0,n):
        for j in range(0,k+1):
            A[k,j] = n-k
        for i in range(0,k):
            A[i,k] = n-k

    resultados = aa.AutovalsAutovecs(A, epsilon, deslocamentos)

    print("Foram necessárias k = {0} iterações no método QR\n".format(resultados[2]))

    for i in range(0, A.shape[0]):
        # Configuração de impressão para ter mais dígitos
        #np.set_printoptions(formatter={'float': '{: 12.10f}'.format})

        # Impressão dos autovals/autovecs
        print('Autovalor: {0:12.10f}, Autovetor:'.format(resultados[0][i]), resultados[1][:, i])
    
    print("\n")
    print('############\n')
    EstimativasErroGerais(A, resultados[0], resultados[1])
    print('############')
    EstimativasErroAnalitico(resultados[0])
    print('############')

teste_b()
