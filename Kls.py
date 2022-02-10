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

def getMasses(decay):
    PDFdir = Path(decay + " PDFs/")
    masses = []
    for x in PDFdir.iterdir():
        if str(x).endswith('.csv'):
            mass = float(str(x)[len(str(PDFdir))+1:-4])
            masses.append(mass)
    return np.sort(np.array(masses))

def archivePDF(x, PDF, DM):
    PDFdir = Path("Kls PDFs/")
    PDFfile = PDFdir / (str(DM) + ".csv")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    with open(PDFfile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(x)
        writer.writerow(PDF)

def archivePlot(x, PDF, DM):
    PDFdir = Path("Kls PDFs/")
    PDFfile = PDFdir / (str(DM) + ".png")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    fig, ax = plt.subplots()
    ax.plot(x, PDF, 'k.')
    ax.set_title('DM(' + str(DM) + ' MeV) -> gamma via Kls')
    ax.set_ylabel('dN/dE (1/Mev)')
    ax.set_xlabel('E (Mev)')
    ax.legend()
    fig.savefig(PDFfile)
    fig.show()

def main():
    masses = getMasses("Kls")
    #trueNorm = 0.7599
    trueNorm = 4*0.307+6*0.195+2*0.125
    for DM in masses:
        x,y = getPDF("Kls", DM)
        '''
        print("origional error in norm for DM "+str(DM)+":")
        print(integrate.simps(y,x)-trueNorm)
        badNorm = integrate.simps(y,x)
        print("badNorm:")
        print(badNorm)
        y = y*trueNorm/badNorm
        print("corrected error in norm for DM "+str(DM)+":")
        print(integrate.simps(y,x)-trueNorm)
        archivePDF(x, y, DM)
        '''
        archivePlot(x, y, DM)

main()
