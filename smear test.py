import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import csv
from pathlib import Path
import os
import random
random.seed()

DM = 1000.0
eta = 548
pi = 135
pipm = 140

'If changing starting energy, CDF, randomE, and background starting energies will change'
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
    return x,y

def smearPDF(energy, decay, DM, res):
    E, PDF = getPDF(decay, DM)
    f = interp1d(E, PDF, kind='linear', fill_value='extrapolate')
    E = np.linspace(startE, max(E), 1000)
    PDF = f(E)/np.sqrt(2*np.pi*res**2)*np.exp(-.5*(energy - E)**2/res**2)
    return integrate.simps(PDF, E)

def getSmearedPDF(decay, DM, res):
    E, PDF = getPDF(decay, DM)
    E2 = np.linspace(startE, max(E) + 4*res, 1000)
    PDF2 = []
    for e in E2:
        PDF2.append(smearPDF(e, decay, DM, res))
    return E2, PDF2

def main():
    DM = 1000.
    x, y = getPDF('Pi Pi Eta', DM)
    a, b = getSmearedPDF('Pi Pi Eta', DM, 30)
    f = interp1d(x, y, kind='linear', fill_value='extrapolate')
    E1 = np.linspace(startE, max(x), 1000)
    fig, ax = plt.subplots()
    ax.plot(E1, f(E1), '-', label='True PDF')
    ax.plot(a, b, '-', label='smeared PDF')
    ax.set_title('1000 MeV Pi Pi Eta deacay spectrum smeared with 30 MeV resolution')
    ax.set_ylabel('counts')
    ax.set_xlabel('Energy [MeV]')
    ax.legend()
    fig.show()

main()
