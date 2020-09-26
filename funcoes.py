# -*- coding: utf-8 -*-

#Importa Bibliotecas Necessárias
from math import log10, pi, sin, cos, ceil, sqrt
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np


class butterworth:
    
    # Essa função inicializa os valores do filtro
    def __init__(self, tipo, Wp1, Ws1, Ap, As, Wp2 = 0, Ws2 = 0):
        self.tipo = tipo
        self.Wp = Wp1
        self.Ws = Ws1
        self.Ap = Ap
        self.As = As
        if tipo == "PF" or tipo == "RF":
            self.Wp1 = Wp1
            self.Wp2 = Wp2
            self.Ws1 = Ws1
            self.Ws2 = Ws2
    
    
    # Bandas de Passagem
    def bandas(self):
        Bp = self.Wp2 - self.Wp1
        Bs = self.Ws2 - self.Ws1
        self.Bs = Bs
        self.Bp = Bp
        return Bp, Bs
    
    
    # Essa função define e retorna a ordem do filtro
    def ordem(self):
        butterworth.bandas(self)
        if self.tipo == "PB":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Ws/self.Wp))
        elif self.tipo == "PA":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Wp/self.Ws))
        elif self.tipo == "PF" or self.tipo == "RF":
            n1 = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Ws2/self.Wp2))
            n2 = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Wp1/self.Ws1))
            n = max(n1, n2)
        N = ceil(n)
        self.N = N
        return n, N
        
    
    # Essa função define e retorna a frequência de corte do filtro
    def freq_corte(self):
        if self.tipo == "PB":
            Wc = self.Wp / pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            self.Wc = Wc
            return Wc
        elif self.tipo == "PA":
            Wc = self.Wp * pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            self.Wc = Wc
            return Wc
        elif self.tipo == "PF":
            Wc1 = self.Wp1 * pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            Wc2 = self.Wp2 / pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            self.Wc1 = Wc1
            self.Wc2 = Wc2
            self.Wc = Wc2
            return Wc1, Wc2
        elif self.tipo == "RF":
            Wc1 = self.Wp1 / pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            Wc2 = self.Wp2 * pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            self.Wc1 = Wc1
            self.Wc2 = Wc2
            self.Wc = Wc1
            return Wc1, Wc2
    
    
    # Frequência de Ressonância
    def freq_ress(self):
        Wo = sqrt(self.Wc1*self.Wc2)
        self.Wo = Wo
        return Wo
    
    
    # Frequência de Ressonância
    def banda_corte(self):
        Bw = self.Wc2 - self.Wc1
        self.Bw = Bw
        return Bw
    
    
    # Essa função define e retorna as raízes do denominador da FT
    def raizes_unit(self):
        Sk = list()
        for k in range(1, self.N + 1):
            Sk.append( complex(-sin( (pi*(2*k-1)) / (2*self.N) ), cos( 
                (pi*(2*k-1)) / (2*self.N) )))
        self.Sk = Sk
        return Sk
    
    
    # Essa função define e retorna a FT do filtro
    def func_tranf(self):
        butterworth.raizes_unit(self)
        butterworth.banda_corte(self)
        poli = list()
        poli = np.poly(self.Sk)
        coefReal = poli.real
        coef = list()
        aux = 0
        for i in range(-self.N, 1):
            coef.append(coefReal[aux]*pow(self.Wc, i))
            aux = aux + 1
        
        if self.tipo == "PB":
            num, den = signal.lp2lp(coefReal[-1], coefReal, self.Wc)
            H = signal.TransferFunction(num, den)
        elif self.tipo == "PA":
            num, den = signal.lp2hp(coefReal[-1], coefReal, self.Wc)
            H = signal.TransferFunction(num, den)
        elif self.tipo == "PF":
            num, den = signal.lp2bp(coefReal[-1], coefReal, self.Wo, self.Bw)
            H = signal.TransferFunction(num, den)
        elif self.tipo == "RF":
            num, den = signal.lp2bs(coefReal[-1], coefReal, self.Wo, self.Bw)
            H = signal.TransferFunction(num, den)
        self.H = H
        return H
    
    
    # Essa função plota o Diagrama de Bode
    def plotar(self):
        # Plotagem do Módulo
        w, y, phase = self.H.bode(w = np.arange(0, 15000, step = 1))
        plt.figure(1)
        plt.grid(True)
        plt.xlim(0, 12000)
        plt.ylim(-70, 0)
        plt.plot(w, y)
    
        # Plotagem da Fase
        plt.figure(2)
        plt.grid(True)
        plt.xlim(0, 6000)
        plt.ylim(-360, 0)
        plt.plot(w, phase, 'r')
            
    
    # Essa função define e retorna os componenetes e seus valores do filtro
    def componentes(self, topologia, Ki):
        elementos = 2*abs(np.real(self.Sk))
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    