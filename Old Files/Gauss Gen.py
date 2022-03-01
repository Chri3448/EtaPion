import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import csv
from pathlib import Path
import os
startE = 10

def archivePlot(x, PDF):
    PDFdir = Path("Gauss4 PDFs/")
    PDFfile = PDFdir / ("1000.0.png")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    fig, ax = plt.subplots()
    ax.plot(x, PDF, 'k.')
    ax.set_title('Gauss4')
    ax.set_ylabel('dN/dE (1/Mev)')
    ax.set_xlabel('E (Mev)')
    ax.legend()
    fig.savefig(PDFfile)

def archivePDF(x, PDF):
    PDFdir = Path("Gauss4 PDFs/")
    PDFfile = PDFdir / ("1000.0.csv")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    with open(PDFfile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(x)
        writer.writerow(PDF)

x= np.linspace(startE, 1000, 1000)
PDF = 1/np.sqrt(2*np.pi*(30)**2)*np.exp(-.5*(x-200)**2/(30)**2)

archivePlot(x, PDF)
archivePDF(x, PDF)
