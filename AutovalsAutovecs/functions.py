import numpy as np
import math as math

def HouseHolder(A, H, n):
    A = np.copy(A)
    H = np.copy(H)
    def dot_prod(a,b,k):
        soma = 0
        for i in range(k,a.shape[0]):
            soma = soma + a[i]*b[i]
        return soma
    
    w = np.zeros(n)
    a = np.zeros(n)

    produto_escalar = np.zeros(n)

    for I in range(0,n-2):
        a = np.array(A[:,I]) # I-ésima coluna de A
        a[I] = 0 # aI = [0, bar_aI]
        w = a+0 # wI = aI
        w[I+1] = w[I+1]+np.sign(a[I+1])*math.sqrt(dot_prod(a,a,I+1)) # aplicando o e da fórmula do w

        for i in range(I+1,n):
            if i>I+1:
                # Aplicando zeros na I-ésima linha e coluna
                A[i,I] = 0
                A[I,i] = 0
            else:
                # Elemento não nulo da I-ésima linha e coluna
                A[i,I] = a[i]-2*dot_prod(w,a,I)/dot_prod(w,w,I)*w[i]
                A[I,i] = a[i]-2*dot_prod(w,a,I)/dot_prod(w,w,I)*w[i]
        # Iterando sobre matriz AI
        # -> Multiplicação por HwI pela esquerda
        for j in range(I+1,n):
            produto_escalar[0] = -2*dot_prod(w,np.squeeze((np.asarray(A[:,j]))),I+1)/dot_prod(w,w,I+1)
            for i in range(I+1,n):
                A[i,j] = A[i,j]+produto_escalar[0]*w[i]
    
        # -> Multiplicação por HwI pela direita
        for j in range(I+1,n):
            produto_escalar[j] = -2*dot_prod(w,np.squeeze((np.asarray(A[j,:]))),I+1)/dot_prod(w,w,I+1)
        for j in range(I+1,n):
            for i in range(j,n):
                A[i,j] = A[j,i]+produto_escalar[j]*w[i]
                A[j,i] = A[i,j]
        # Multiplicação para encontrar HT
        for i in range(0,n):
            produto_escalar[i] = -2*dot_prod(w,np.squeeze((np.asarray(H[:,i]))),I+1)/dot_prod(w,w,I+1)
        for i in range(0,n):
            for j in range(I+1,n):
                H[i,j]=H[i,j]+produto_escalar[i]*w[j]
    return (A, H)

def metodo_QR(alfa, beta, epsilon, deslocamentos, V):
    # Este método faz o cálculo dos cossenos (forma estável)
    def __calcula_c(alfa, beta, i):
        if abs(alfa[i]) > abs(beta[i]):
            tau = -beta[i] / alfa[i]
            return 1 / math.sqrt(1 + tau * tau)
        else:
            tau = -alfa[i] / beta[i]
            return tau / math.sqrt(1 + tau * tau)

    # Este método faz o cálculo dos senos (forma estável)
    def __calcula_s(alfa, beta, i):
        if abs(alfa[i]) > abs(beta[i]):
            tau = -beta[i] / alfa[i]
            return tau / math.sqrt(1 + tau * tau)
        else:
            tau = -alfa[i] / beta[i]
            return 1 / math.sqrt(1 + tau * tau)

    # Este método faz o cálculo de mu
    # m: o tamanho-1 da submatriz onde o método é aplicado
    def __calcula_mu(alfa, beta, m, deslocamentos):
        # O if abaixo detecta se foram solicitados deslocamentos espectrais
        if deslocamentos:
            # Precisamos obter primeiramente dk.
            dk = (alfa[m - 1] - alfa[m]) / 2
            mu = alfa[m] + dk - np.sign(dk) * math.sqrt(dk ** 2 + beta[m - 1] ** 2)
            return mu
        else:
            return 0

    # Este método faz a normalização de vector
    def __normaliza(vector):
        soma = 0
        for i in range(0, len(vector)):
            soma = soma + vector[i] ** 2
        return vector / math.sqrt(soma)

    # Este é o método público desta classe
    c = np.zeros(alfa.shape[0] - 1) # Vetor que armazena os valores de c em cada subiteração i da decomposição
                                    # QR da k-ésima iteração
    s = np.zeros(alfa.shape[0] - 1) # Vetor que armazena os valores de s em cada subiteração i da decomposição
                                    # QR da k-ésima iteração
    mu = 0

    # De acordo com o for, esses 3 vetores mudam de função
    Diagonal_principal = np.copy(alfa)
    Subdiagonal_superior = np.copy(beta)
    Subdiagonal_inferior = np.copy(beta)

    # Este vetor é usado para armazenar o valor intermediário de V para que ele seja atribuído no final
    temp_V = np.copy(V)

    # O algoritmo itera sobre matriz m+1 por m+1, que vai diminuindo à medida que os autovalores são encontrados
    m = alfa.shape[0] - 1

    k = 0 # Contador de iterações

    while (m >= 0):
        # Estes quatro vetores temp são usados como temporários para armazenar as células da
        # matriz que são atualizadas entre subiterações dos três fors
        temp_diagonal_principal = np.zeros(2)
        temp_subdiagonal = np.zeros(2)

            # Quando multiplicamos o subresultado de V(k+1) por Qi(k), só são modificadas
            # duas colunas do subresultado, aqui estão os dois temporários das duas colunas modificadas
        temp_V_coluna_esquerda = np.zeros(alfa.shape[0])
        temp_V_coluna_direita = np.zeros(alfa.shape[0])

        # O if só é executado se a submatriz é pelo menos 2x2 e o último beta não puder ser zerado
        if (m > 0 and abs(Subdiagonal_inferior[m - 1]) >= epsilon):
            if (k > 0):
                mu = __calcula_mu(Diagonal_principal, Subdiagonal_inferior, m, deslocamentos)
                Diagonal_principal = Diagonal_principal - mu
            # Este for faz a decomposição QR
            #
            # Este for realiza m subiterações i
            # A cada subiteração o subresultado de R(k) é multiplicado por Qi(k)
            # Até no final ter a própria R(k)
            #
            # Funções dos seguinte vetores neste for:
            # Diagonal_principal: Vetor que armazena a diagonal principal do subresultado de R(k),
            # e. no final do for, da própria R
            #
            # Subdiagonal_superior: Vetor que armazena a 1ª subdiagonal acima da diagonal principal do
            # subresultado de R(k), e. no final do for, da própria R(k)
            #
            # Subdiagonal_inferior: Vetor que armazena a 1ª subdiagonal abaixo da diagonal principal do
            # subresultado de R(k), e. no final do for, da própria R(k)
            #
            #
            #
            # Otimizações implementadas neste for:
            # 1. Só são armazenados os cossenos e senos, e não a matriz Qi(k) inteira
            # 2. O elemento da subiteração i do subresultado de R(k) de posição (i,i+1) é zerado automaticamente
            # 3. Só se calcula os valores da diagonal principal e das duas subdiagonais, imediatamente abaixo e
            # acima da principal, pois os elementos das outras diagonais são zerados e não são usados
            # no cálculo de A(k+1)
            # 4. Entre cada subiteração somente 5 células importantes são alteradas no subresultado de R(k),
            # exceto na última iteração
            for i in range(0, m):
                # Calcula c e s da subiteração i da decomposição QR #Otimização 1
                c[i] = __calcula_c(Diagonal_principal, Subdiagonal_inferior, i)
                s[i] = __calcula_s(Diagonal_principal, Subdiagonal_inferior, i)
                # Calcula os novos valores das duas células atualizadas da diagonal principal
                # do subresultado de R(k)
                temp_diagonal_principal[0] = c[i] * Diagonal_principal[i] - s[i] * Subdiagonal_inferior[i]
                temp_diagonal_principal[1] = s[i] * Subdiagonal_superior[i] + c[i] * Diagonal_principal[i + 1]

                # Calcula os novos valores das duas células atualizadas da 1ª subdiagonal acima
                # da diagonal principal do subresultado de R(k)
                temp_subdiagonal[0] = c[i] * Subdiagonal_superior[i] - s[i] * Diagonal_principal[i + 1]
                if (i != m - 1):
                    # Não existe Subdiagonal_superior[i+1] na última subiteração, por isso este if
                    temp_subdiagonal[1] = c[i] * Subdiagonal_superior[i + 1]

                # Atualiza o subresultado com os valores temporários
                Subdiagonal_inferior[i] = 0.0  # Otimização 2
                Diagonal_principal[i] = np.copy(temp_diagonal_principal[0])
                Diagonal_principal[i + 1] = np.copy(temp_diagonal_principal[1])
                Subdiagonal_superior[i] = np.copy(temp_subdiagonal[0])

                if (i != m - 1):
                    Subdiagonal_superior[i + 1] = np.copy(temp_subdiagonal[1])
            
            # Este for faz o cálculo de A(k+1)
            #
            # Este for realiza m subiterações i
            # A cada subiteração o subresultado de A(k+1) é multiplicado por Qi(k)T
            # Até no final ter a própria A(k+1)
            #
            # Neste for
            # Diagonal_principal: Vetor que armazena a diagonal principal do subresultado de A(k+1),
            # e. no final do for, da própria A(k+1)
            #
            # Subdiagonal_superior: Vetor que armazena a 1ª subdiagonal acima da diagonal principal
            # do subresultado de A(k+1), e. no final do for, da própria A(k+1)
            #
            # Subdiagonal_inferior: Vetor que armazena a 1ª subdiagonal abaixo da diagonal principal
            # do subresultado de A(k+1), e. no final do for, da própria A(k+1)
            #
            #
            #
            # Otimizações implementadas neste for:
            # 1. Como a matriz A(k+1) termina simétrica, então só calculamos o valor para a
            # diagonal imediatamente abaixo da principal
            # 2. Como as multiplicações por Qi(k)T só modificam duas colunas do subresultado de A(k+1),
            # então só é necessário calcular 3 células do subresultado
            for i in range(0, m):
                # Calcula os novos valores das duas células atualizadas
                # da diagonal principal do subresultado de A(k+1)
                temp_diagonal_principal[0] = c[i] * Diagonal_principal[i] - s[i] * Subdiagonal_superior[i]
                temp_diagonal_principal[1] = c[i] * Diagonal_principal[i + 1]

                # Calcula o novo valos da célula atualizada da 1ª subdiagonal abaixo
                # da diagonal principal do subresultado de A(k+1)
                temp_subdiagonal[0] = -s[i] * Diagonal_principal[i + 1]  # Otimização 1

                # Atualiza o subresultado com os valores temporários
                Diagonal_principal[i] = np.copy(temp_diagonal_principal[0])
                Diagonal_principal[i + 1] = np.copy(temp_diagonal_principal[1])
                Subdiagonal_inferior[i] = np.copy(temp_subdiagonal[0])

                # Como A(k+1) é simétrica, não precisamos calcular a subdiagonal superior #Otimização 1
                Subdiagonal_superior[i] = np.copy(temp_subdiagonal[0])

            # Faz A(k+1) = R(k)*Q(k)+muk*I
            if (k > 0):
                Diagonal_principal = Diagonal_principal + mu

            # Este for faz a multiplicação V(k+1)=V(k)*Q(k)T, e armazena o resultado em temp_V
            #
            # Este for realiza m subiterações subiteracao
            # A cada subiteração o subresultado de V(k+1) é multiplicado por Qi(k)T
            # Até no final ter a própria V(k+1)
            #
            # Otimizações implementadas neste for:
            # 1. Como as multiplicações por Qi(k)T só modificam duas colunas do subresultado de V(k+1),
            # então só é necessário calcular 2 colunas do subresultado
            for subiteracao in range(0, m):
                for linha in range(0, alfa.shape[0]):
                    temp_V_coluna_esquerda[linha] = temp_V[linha, subiteracao] * c[subiteracao]\
                                                    - temp_V[linha, subiteracao + 1] * s[subiteracao]
                    temp_V_coluna_direita[linha] = temp_V[linha, subiteracao + 1] * c[subiteracao]\
                                                    + temp_V[linha, subiteracao] * s[subiteracao]

                    temp_V[linha, subiteracao] = temp_V_coluna_esquerda[linha]
                    temp_V[linha, subiteracao + 1] = temp_V_coluna_direita[linha]

            k = k + 1


        # Este if faz a eliminação de beta se este for menor que epsilon e diminui o escopo (m=m-1)
        # do algoritmo, para que ele passe a trabalhar com a submatriz
        if (abs(Subdiagonal_inferior[m - 1]) < epsilon):
            Subdiagonal_inferior[m - 1] = 0.0
            Subdiagonal_superior[m - 1] = 0.0
            m = m - 1

    # __normaliza autovetores
    for j in range(0, alfa.shape[0]):
        temp_V[:, j] = __normaliza(temp_V[:, j])
    return (Diagonal_principal, temp_V, k)

def AutovalsAutovecs(A, epsilon, deslocamentos):
    H = np.identity(4)
    householder = HouseHolder(A, H, A.shape[0])
    alfa = np.zeros(A.shape[0])
    beta = np.zeros(A.shape[0]-1)

    for i in range(0, householder[0].shape[0]):
        alfa[i] = householder[0][i,i]
    for i in range(0, householder[0].shape[0]-1):
        beta[i] = householder[0][i+1,i]
    
    return metodo_QR(alfa, beta, epsilon, deslocamentos, householder[1])


