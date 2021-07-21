import AutovalsAutovecs as aa
import numpy as np
import math as math

def EstimativasErroGerais(A, Lambda, HV):
    n = A.shape[0]

    print("===Estimativas de Erro Gerais===\n")

    # Faz uma estimativa de erro = max(abs(H*V*VT*HT-I))
    erro = np.max(np.abs(np.matmul(HV, np.transpose(HV)) - np.identity(n)))
    print("(Verificação de ortogonalidade)")
    print('max(abs((HT*V)*(HT*V)T-I)) = {0}'.format(erro))

    # Faz uma estimativa de erro = max(abs(matmul(A,H*V)-matmul(H*V,L)))
    L = np.zeros((n, n))

    for i in range(0, n):
        L[i, i] = Lambda[i]

    erro = np.max(np.abs(np.matmul(A, HV) - np.matmul(HV, L)))
    print('\nmax(abs(matmul(A,HT*V)-matmul(HT*V,L))) = {0}\n'.format(erro))

    # Faz uma estimativa de erro comparando (A*v)/v, com os autovalores obtidos
    erros = np.zeros(n)
    for i in range(0, n):
        erros[i] = math.sqrt(math.pow(np.mean(np.divide(np.matmul(A, HV[:, i]), HV[:, i])) - Lambda[i], 2))
        print('Autovalor obtido fazendo matdiv(matmul(A,v),v): {0:12.10f},'
                '  Autovalor obtido pelo método: {1:12.10f}, erro: {2}'
                .format(np.mean(np.divide(np.matmul(A, HV[:, i]), HV[:, i])), Lambda[i], erros[i]))
    print('Erro máximo = {0}\n'.format(np.max(erros)))

    # Faz estimativa de erro verificando o quão próximo de zero está A*v-lambda*v
    print("\nEstimativa de erro fazendo max(abs(A*v-lambda*v)) para cada autovalor")
    for i in range(0, n):
        erros[i] = np.max(np.abs(A.dot(HV[:, i])-Lambda[i]*HV[:, i]))
        print('Autovalor: {0:12.10f}, Erro {1:12.10f}'
                .format(Lambda[i], erros[i]))
    print('Erro máximo = {0}\n'.format(np.max(erros)))

def teste_a():
    # Menu de seleção
    print('\n############')
    print('Teste a selecionado')
    print('############\n')

    epsilon = 0.1
    deslocamentos = True

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
    
    print("\n======================RESULTADOS DO TESTE A======================")

    A = np.matrix([[2,4,1,1],[4,2,1,1],[1,1,1,2],[1,1,2,1.0]])

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

