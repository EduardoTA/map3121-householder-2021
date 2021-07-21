import AutovalsAutovecs as aa
import numpy as np
np.set_printoptions(linewidth=500)
np.set_printoptions(formatter={'float': '{: 12.10f}'.format})
np.set_printoptions(precision=10)
import math as math

# def check_symmetric(a, tol=1e-8):
#     return np.all(np.abs(a-a.T) < tol)

def teste_c():
    print('\n############')
    print('Teste c selecionado')
    print('############\n')

    epsilon = 0.1
    deslocamentos = True
    f = 0 # Arquivo será inserido nessa variável

    n_total = 0 # Número de nós da treliça
    n_moveis = 0 # Número de nós móveis da treliça
    n_barras = 0 # Número de barras da treliça

    rho = 0 # Densidade rho [kg/m³]
    A = 0 # Área transversal das barras [m²]
    E = 0 # Módulo de elasticidade em Pa

    K = 0 # Matrix de rigidez total
    m = 0 # Vetor de massas dos nós
    M = 0 # Matriz onde as massas dos nós estarão na diagonal principal

    # Fazendo leitura das entradas de usuário
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

    ## Faz a leitura das duas primeiras linhas do arquivo
    try:
        n_total = int(f.pop(0)) # Número de nós da treliça
        n_moveis = int(f.pop(0)) # Número de nós móveis da treliça
        n_barras = int(f.pop(0)) # Número de barras da treliça

        rho = float(f.pop(0)) # Densidade rho [kg/m³]
        A = float(f.pop(0)) # Área transversal das barras [m²]
        E = float(f.pop(0)) # Módulo de elasticidade em GPa
        E = E*1e9 # Módulo de elasticidade em Pa

        K = np.zeros((n_moveis*2, n_moveis*2)) # Matrix de rigidez total
        m = np.zeros(n_moveis) # Vetor de massas dos nós
        M = np.zeros((n_moveis*2, n_moveis*2)) # Matriz onde as massas dos nós estarão na diagonal principal
    
    except:
        print('Erro na leitura dos parâmetros das 2 primeiras linhas do arquivo')
        return
        

    # Uma iteração é realizada para cada barra
    for barra in range(0, n_barras):
        try:
            # Extraindo cada linha do arquivo correspondente às barras
            i = int(f[barra*4]) # Extrai o nó i
            j = int(f[barra*4+1]) # Extrai o nó j
            if i > n_total or j > n_total: # Se for detectado i ou j inválido, levantar exceção
                raise
            i = i-1 # Faz o índice começar em 0
            j = j-1 # Faz o índice começar em 0
            theta = float(f[barra*4+2]) # Extrai o ângulo da barra
            L = float(f[barra*4+3]) # Extrai o comprimento da barra
        except:
            print('Erro na leitura dos parâmetros das barras da treliça')
            return

        alfa = A*E/L # Coeficiente da equação
        C = math.cos(theta*math.pi/180) # Cosseno da equação
        S = math.sin(theta*math.pi/180) # Seno da equação

        # Fazendo o cálculo de cada elemento da matriz da equação
        # E somando à matriz K na posição correta
        n = n_moveis*2 # Dimensão das matrizes K e M
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

        # Nestes dois ifs, as massas de cada nó são acumuladas no vetor m
        if 0 <= i < n_moveis:
            m[i] = m[i] + 1/2*rho*A*L
        if 0 <= j < n_moveis:
            m[j] = m[j] + 1/2*rho*A*L

    
    # Agora a matriz K e o vetor m estão completamente montados

    # Agora serão montadas as matrizes M, tilde_M e tilde_K

    tilde_M = np.zeros((n_moveis*2, n_moveis*2)) # Instancia tilde_M = M^(-1/2)
    
    # Este for monta as matrizes M e tilde_M
    for no in range(0, n_moveis):
        M[no*2,no*2] = m[no]
        M[no*2+1, no*2+1] = m[no]

        tilde_M[no*2,no*2] = 1/math.sqrt(m[no])
        tilde_M[no*2+1, no*2+1] = tilde_M[no*2,no*2]

    # Esta multiplicação de matrizes monta a matriz tilde_K
    tilde_K = np.matmul(np.matmul(tilde_M, K), tilde_M) # M^(-1/2)*K*M^(-1/2)

    # Aqui o método QR com Householder é aplicado
    resultados = aa.AutovalsAutovecs(tilde_K, epsilon, deslocamentos)

    # Calcula as frequências
    w = np.zeros(n_moveis*2)
    for i in range(n_moveis*2):
        w[i] = math.sqrt(resultados[0][i]) # w = sqrt(lambda)

    # Calcula os modos de vibração
    z = np.zeros((n_moveis*2, n_moveis*2))
    for j in range(n_moveis*2):
        z[:,j] = tilde_M.dot(resultados[1][:,j]) # z = M^(-1/2)*y
    
    # Monta uma lista com as frequencias e modos naturais
    list = []
    for i in range(0, n_moveis*2):
        list.append([w[i], z[:, i]])
    
    list = sorted(list,key=lambda x: x[0]) # Ordena essa lista da menor freq. para a maior

    print("\n======================RESULTADOS DO TESTE C======================")
    
    print("\nForam necessárias k = {0} iterações do método QR\n".format(resultados[2]))

    # Seletor de opções:
    def seletor():
        print('Selecione uma opção:')
        print('(0) Para imprimir as frequências e modos naturais no terminal')
        print('(1) Dump de matrizes e vetores do programa em arquivos')
        print('(2) Estimativas de erro')
        print('(qualquer outra) Para sair do teste')
        selecao = input('= ')

        if selecao == '0':
            # Esta seleção faz a apenas a impressão no terminal das freqs. e modos naturais
            print('\n')
            for i in range(0, n_moveis*2):
                # Configuração de impressão para ter mais dígitos
                np.set_printoptions(formatter={'float': '{: 12.10f}'.format})

                # Impressão dos autovals/autovecs
                print('freq: {0:12.10f},\nmodo:'.format(list[i][0]), list[i][1])
                print('\n')
            seletor()
        elif selecao == '1':
            # Esta seleção imprime as matrizes em arquivos
            print('\n')
            try:
                np.savetxt('m.txt', m)
                np.savetxt('K.txt', K)
                np.savetxt('tilde_K.txt', tilde_K)
                np.savetxt('M.txt', M)
                np.savetxt('tilde_M.txt', tilde_M)
                np.savetxt('w.txt', w)
                np.savetxt('z.txt', z)
                print('Arquivos criados com sucesso')
            except:
                print('Arquivos não foram criados com sucesso')
            seletor()
        elif selecao == '2':
            # Esta seleção faz estimativas de erro
            print('\n')
            print('Estimativa de erro comparando sqrt((K*z)/(M*z)) com as frequências obtidas')
            erros = np.zeros(n)
            for i in range(0, n_moveis*2):
                erros[i] = np.abs(np.sqrt(np.max(np.divide(K.dot(list[i][1]), M.dot(list[i][1])))) - list[i][0])
                print('Freq. obtida fazendo sqrt((K*z)/(M*z)): {0:12.10f},'
                        '  Freq. obtida pelo método: {1:12.10f}, erro: {2}'
                        .format(np.sqrt(np.max(np.divide(K.dot(list[i][1]), M.dot(list[i][1])))), list[i][0], erros[i]))
            print('Erro máximo = {0}\n'.format(np.max(erros)))
            
            print('\n')
            print('Estimativa de erro calculando M*x\" + K*x, para t=0s')

            erros = np.zeros(n)
            for i in range(0, n_moveis*2):
                erros[i] = np.max(np.abs(-(list[i][0]**2)*M.dot(list[i][1])+K.dot(list[i][1])))
                print('freq.: {0:12.10f}, |-w^2*M*z + K*z|: {1:12.10f}'
                        .format(list[i][0],
                        np.max(np.abs(-(list[i][0]**2)*M.dot(list[i][1])+K.dot(list[i][1]))))
                    )
            print('Erro máximo = {0}\n'.format(np.max(erros)))
            print('\n')
            seletor()
        else:
            return
    seletor()
    return
