# -*- coding: utf-8 -*-

#Importa Bibliotecas Necessárias
from math import log10, pi, sin, cos, ceil
from scipy import signal
import matplotlib.pyplot as plt
import numpy


class butterworth:
    
    # Essa função inicializa os valores do filtro
    def __init__(self, tipo, Wp, Ws, Ap, As):
        self.tipo = tipo
        self.Wp = Wp
        self.Ws = Ws
        self.Ap = Ap
        self.As = As
    
    # Essa função define e retorna a ordem do filtro
    def ordem(self):
        
        if self.tipo == "PB":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Ws/self.Wp))
            N = ceil(n)
        elif self.tipo == "PA":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Wp/self.Ws))
            N = ceil(n)
        self.N = N
        return n, N
        
    def freq_corte(self):
        if self.tipo == "PB":
            Wc = self.Wp / pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
        elif self.tipo == "PA":
            Wc = self.Wp * pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
        self.Wc = Wc
        return Wc
    
    def raizes_unit(self):
        Sk = list()
        for k in range(1, self.N + 1):
            Sk.append( complex(-sin( (pi*(2*k-1)) / (2*self.N) ), cos( 
                (pi*(2*k-1)) / (2*self.N) )))
        self.Sk = Sk
        return Sk
    
    def func_tranf(self):
        butterworth.raizes_unit(self)
        poli = list()
        poli = numpy.poly(self.Sk)
        coefReal = poli.real
        coef = list()
        if self.tipo == "PB":
            aux = 0
            for i in range(-self.N, 1):
                coef.append(coefReal[aux]*pow(self.Wc, i))
                aux = aux + 1
            H = signal.TransferFunction(coef[-1], coef)
        elif self.tipo == "PA":
            den = list()
            aux = self.N
            for i in range(0, self.N+1):
                coef.append(coefReal[aux]*pow(self.Wc, i))
                aux = aux - 1
                if i == 0:
                    den.append(coef[0])
                else:
                    den.append(0)
            H = signal.TransferFunction(den, coef)
        self.H = H
        return H
    
    def plotar(self):
        # Plotagem do Módulo
        w, y, phase = signal.bode(self.H)
        plt.figure(1)
        plt.grid(True)
        plt.plot(w, y)
        
        # Plotagem da Fase
        plt.figure(2)
        plt.grid(True)
        plt.plot(w, phase, 'r')
            
    def componentes(self, topologia, Ki):
        elementos = 2*abs(numpy.real(self.Sk))
        comp = list()
        
        # Definição da topologia utilizada
        if topologia == "LCL":
            for i in range(0, self.N):
                if i%2 == 0:
                    comp.append("L")
                else:
                    comp.append("C")
        elif topologia == "CLC":
            for i in range(0, self.N):
                if i%2 == 0:
                    comp.append("C")
                else:
                    comp.append("L")
        self.comp = comp
        
        # Transformação de Frequência
        # Para um Passa Baixa
        if self.tipo == "PB":
            elementos = elementos / self.Wc
        # Para um Passa Alta
        elif self.tipo == "PA":
            elementos = 1 / (elementos * self.Wc)
        
        # Transformação de Impedância
        R = Ki
        aux = 0
        for i in comp:
            if i == "C":
                elementos[aux] = elementos[aux] / Ki
            elif i == "L":
                elementos[aux] = elementos[aux] * Ki
            aux = aux + 1
        
        # Definição das Variáveis e Retorno
        self.R = R
        self.elementos = elementos
        return elementos, comp, R
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    