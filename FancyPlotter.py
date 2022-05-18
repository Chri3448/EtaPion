import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import scipy.optimize
import scipy.stats
import scipy.special
import csv
from pathlib import Path
import os
import random
random.seed()
import time
import pickle as pk

startE = 10

def getPDF(decay, DM):
    x = []
    y = []
    PDFdir = Path(decay + " PDFs/")
    PDFfile = PDFdir / (str(DM) + ".csv")
    with open(PDFfile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        for row in reader:
            if x == []:
                for val in row:
                    x.append(float(val))
                continue
            if y == []:
                for val in row:
                    y.append(float(val))
    return np.array(x), np.array(y)

def smearPDF(energy, E, PDF, res):
    sig = res/(2*np.sqrt(2*np.log(2)))
    f = interp1d(E, PDF, kind='linear', fill_value='extrapolate')
    E = np.linspace(min(E), max(E), 1000)
    E_values = np.ones((np.size(energy), np.size(E), 2))
    E_values[:,:,0] = np.tile(energy, (np.size(E),1)).T
    E_values[:,:,1] = np.tile(E, (np.size(energy),1))
    PDF = f(E_values[:,:,1])/np.sqrt(2*np.pi*(sig*E_values[:,:,1])**2)*np.exp(
        -.5*(E_values[:,:,0] - E_values[:,:,1])**2/(sig*E_values[:,:,1])**2)
    #return integrate.trapz(PDF, E, axis = 1)
    return np.sum(PDF[:,:-1], axis = 1)*(E[1]-E[0])

def getSmearedPDF(decay, DM, res):
    sig = res/(2*np.sqrt(2*np.log(2)))
    E, PDF = getPDF(decay, DM)
    E2 = np.linspace(startE, max(E) + 4*sig*max(E), 1000)
    PDF2 = smearPDF(E2, E, PDF, res)
    return E2, PDF2


def SpectrumPlot(decay, DM):
    x,y = getPDF(decay, DM)
    x = x[np.where(x < 350.0)]
    y = y[0: np.size(x)]
    y = y/integrate.simps(y,x)
    spec = interp1d(x,y, kind='linear')
    x=np.linspace(0.1,np.max(x), 10000)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=21)
    fig, ax = plt.subplots(figsize=(10,7))
    ax.plot(x, spec(x), 'k')
    plt.xlabel(r'$E_\gamma \hspace{.1cm} (MeV)$')
    plt.ylabel(r'$\frac{dN_\gamma}{dE_\gamma} \hspace{.1cm} (MeV^{-1})$')
    if decay == 'Kls':
        plt.title(r'$K_LK_S$')
        #plt.title(r'$K_LK_S--\sqrt{s} = '+str(DM)+' MeV$')
    elif decay == 'K+-':
        plt.title(r'$K^+K^-$')
    elif decay == 'Pi Pi Eta':
        plt.title(r'$\pi\pi\eta$')
    elif decay == 'Pi Eta':
        plt.title(r'$\pi^0\eta$')
    plt.show()
    

def main():
    decay = 'Pi Pi Eta'
    DM = 835.0
    SpectrumPlot(decay, DM)
    
main()
