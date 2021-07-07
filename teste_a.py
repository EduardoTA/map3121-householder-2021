import AutovalsAutovecs as aa
import numpy as np

def teste_a():
    # Menu de seleção
    print('############')
    print('Teste a selecionado')
    print('############')

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
    
    print("\nOs autovalores e autovetores da matriz do item 4.1.a serão impressos\n")

    A = np.matrix([[2,4,1,1],[4,2,1,1],[1,1,1,2],[1,1,2,1.0]])

    resultados = aa.AutovalsAutovecs(A, epsilon, True)

    print("Foram necessárias k = {0} iterações no método QR\n".format(resultados[2]))

    for i in range(0, A.shape[0]):
        # Configuração de impressão para ter mais dígitos
        #np.set_printoptions(formatter={'float': '{: 12.10f}'.format})

        # Impressão dos autovals/autovecs
        print('Autovalor: {0:12.10f}, Autovetor:'.format(resultados[0][i]), resultados[1][:, i])
    
    Lambda = np.zeros((A.shape[0], A.shape[0]))
    for i in range(0, A.shape[0]):
        Lambda[i,i] = resultados[0][i]
    
    erro = np.max(np.abs((np.matmul(np.matmul(resultados[1], Lambda), np.transpose(resultados[1])) - A)))
    print("\nEstimativa de erro H*L*HT-A = {0}\n".format(erro))
    print('############')


teste_a()