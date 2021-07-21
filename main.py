# Nome: Eduardo Tadashi Asato,  nusp: 10823810

from teste_a import teste_a as teste_a
from teste_b import teste_b as teste_b
from teste_c import teste_c as teste_c
from teste_d import teste_d as teste_d

def seletor():
    print('\n============================================')
    print("\nExercício-programa 2 de MAP3121, 2021")
    print("Nome: Eduardo Tadashi Asato,     nusp: 10823810\n")

    print('============================================\n')
    print('Selecione uma opção:')
    print('(0 ou a) Para iniciar o teste a (Cálculo dos autovals e autovecs da matriz no input-a)')
    print('(1 ou b) Para iniciar o teste b (Cálculo dos autovals e autovecs da matriz no input-b)')
    print('(2 ou c) Para iniciar o teste c (Cálculo da treliça)')
    print('(3 ou d) Para iniciar o teste d (Cálculo dos autovals e autovecs de matriz em qualquer arquivo)')
    print('(qualquer outra) Para sair')
    selecao = input('= ')

    if selecao == '0' or selecao == 'a':
        teste_a()
        seletor()
    elif selecao == '1' or selecao == 'b':
        teste_b()
        seletor()
    elif selecao == '2' or selecao == 'c':
        teste_c()
        seletor()
    elif selecao == '3' or selecao == 'd':
        teste_d()
        seletor()
    else:
        return

seletor()
