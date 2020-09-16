from math import log10, pi, sin, cos, ceil
from scipy import signal
import matplotlib.pyplot as plt
import numpy

# Definição de Variáveis
Ap = -1  # dB
As = -40  # dB
Wp = 1000  # rad/s
Ws = 5000  # rad/s
Rl = 50  # Ohms
Rs = Rl  # Ohms // Rl = Rs -> Consideração minha

# Definição da Ordem do Filtro
n = log10( (pow(10, (-As/10)) - 1) / (pow(10,(-Ap/10)) - 1) ) / (2*log10(Ws/Wp))
N = ceil(n)

# Definição da Frequência de Corte
Wc = Wp / pow((pow(10, (-Ap/10))-1), (1/(2*N)))

# Definição dos polos unitários
Sk = list()
for k in range(1, N+1):
    Sk.append( complex(-sin( (pi*(2*k-1)) / (2*N) ), cos( (pi*(2*k-1)) / (2*N) )))
    
# Extrai as raízes da função
raizes = list()
raizes = numpy.poly(Sk)

# Transforma Sk em vetor de coeficientes reais
coefMod = raizes.real

# Definição dos Coeficientes do Denominador
aux = 0
coef = list()
for i in range(-N, 1):
    coef.append(coefMod[aux]*pow(Wc, i))
    aux = aux + 1

# Função de Transferência 
H = signal.TransferFunction(coef[-1], coef)

# Plotagem do Módulo
w, y, phase = signal.bode(H)
plt.figure(1)
plt.grid(True)
plt.plot(w, y)

# Plotagem da Fase
plt.figure(2)
plt.grid(True)
plt.plot(w, phase, 'r')





