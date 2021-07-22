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

#takes as input PDF and returns CDF evaluated at E
def CDF(PDF, E):
    x = np.linspace(startE , E, 100)
    values = []
    for energy in x:
        if not PDF(energy) < 0:
            values.append(PDF(energy))
        else:
            values.append(0)
    return integrate.simps(values, x)

#takes PDF and returns random E(MeV)
def randomE(PDF, norm):
    tolerance = 0.0000000001
    CDFvalue = random.random()*norm
    energyRange = 10
    E = startE + 10
    CDFE = CDF(PDF, E)
    while CDFE < CDFvalue:
        try:
            E +=  10
            energyRange += 10
            CDFE = CDF(PDF, E)
        except:
            E -= 10
            energyRange -= 10
            newE = E
            step = 10
            while newE == E:
                try:
                    CDF(PDF, E + step)
                    newE = E + step
                except:
                    step /= 2
            CDFE = CDF(PDF, newE)
            print(CDFE, CDFvalue)
            E =  newE
            energyRange = newE
    dif = CDF(PDF, E) - CDFvalue
    while np.abs(dif) > tolerance:
        energyRange /= 2
        if dif > 0:
            E -= energyRange
        elif dif < 0:
            E += energyRange
        else:
            break
        dif = CDF(PDF, E) - CDFvalue
    return [E, dif]

#takes PDF and returns statistical ensemble of N particles as list      
def randomEnsemble(PDF, norm, N):
    ensemble = []
    for n in range(N):
        print(n)
        ensemble.append(randomE(PDF, norm)[0])
    return ensemble

#gaussian function used in energy resolution smear
#def gauss(x, mu, sig):
#    return (1/(sig*np.sqrt(2*np.pi)))*np.exp(-.5*(x-mu)**2/sig**2)
    

#taes photon ensemble and energy resolution, returns ensemble smeared by energy resolution
def smearedEnsemble(ensemble, res):
    newEnsemble = []
    for E in ensemble:
        print(E)
        newE = np.random.normal(E, E*res)
        if newE >= startE:
            newEnsemble.append(newE)
    return newEnsemble
        
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

#background PDF as a function of energy
#energy (MeV), exposure (cm^2 yr), Sangle (sr) 
def background(energy, exposure, Sangle):
    secPerYear = 3.154*10**7
    preFactor = 0.00274
    back = preFactor*exposure*Sangle*secPerYear/energy**2
    return back

def backgroundEnsembleGen(exposure, Sangle):
    E = np.linspace(startE, 1000, 10000)
    PDF = background(E, exposure, Sangle)
    specNorm = integrate.simps(PDF, E)
    expectedEvents = int(round(specNorm))
    print(expectedEvents)
    N_Events = np.random.poisson(expectedEvents)
    print(N_Events)
    f = interp1d(E, PDF, kind='linear', fill_value='extrapolate')
    Ensemble = randomEnsemble(f, specNorm, N_Events)
    #sEnsemble = smearedEnsemble(Ensemble, .3)
    file = Path() / ("background.csv")
    with open(file, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(Ensemble)
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(Ensemble, 100, density=True, label='Photon Histogram')
    ax.plot(bins, f(bins)/specNorm, '-', label='PDF')
    ax.legend()
    ax.set_title(str(N_Events) + ' Simulated Photons from background')
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlabel('Energy (Mev)')
    fig.show()

#makes a list of photon energies and saves as .csv
def archiveEnsemblePlot(decay, DM, E, PDF, photons):
    Edir = Path(decay + " Ensembles/")
    Efile = Edir / (str(DM) + ".csv")
    Eplotfile = Edir / (str(DM) + ".png")
    try:
        os.mkdir(Edir)
    except:
        print(Edir,'already exists')
    f = interp1d(E, PDF, kind='linear', fill_value='extrapolate')
    newE = np.linspace(startE, max(E), 1000)
    newPDF = []
    for energy in newE:
        newPDF.append(f(energy))
    normalization = integrate.simps(newPDF, newE)
    Ensemble = randomEnsemble(f, normalization, photons)
    #sEnsemble = smearedEnsemble(Ensemble, 30)
    #print('Photons lost by smear: {}'.format(len(Ensemble)-len(sEnsemble)))
    with open(Efile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(Ensemble)
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(Ensemble, 100, density=True, label='Photon Histogram')
    ax.plot(bins, f(bins)/normalization, '-', label='PDF')
    ax.legend()
    ax.set_title(str(photons) + ' Simulated Photons from DM(' + str(DM) + ' MeV) -> gamma via ' + decay)
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlabel('Energy (Mev)')
    fig.savefig(Eplotfile)
    fig.show()


def main():
    exposure = 30000
    
    '''
    backgroundEnsembleGen(exposure, 0.000404322)
    '''
    

    dRate = 10**-26
    secPerYear = 3.154*10**7
    JfactorDracoDecay = 5.77*(10**24) #(MeV cm^-2 sr^-1)
    
    decay = 'Pi Eta'
    x,y = getPDF(decay, DM)
    f = interp1d(x, y, kind='linear', fill_value='extrapolate')
    E = np.linspace(startE, max(x), 1000)
    PDF = []
    for energy in E:
        PDF.append(f(energy))
    specNorm = integrate.simps(PDF, E)
    expectedEvents = int(round((dRate)/(4*np.pi*DM)*JfactorDracoDecay*specNorm*exposure*secPerYear*0.000404322))
    print('expected events: ' + str(expectedEvents))
    N_Events = np.random.poisson(expectedEvents)
    print(N_Events)

    archiveEnsemblePlot(decay, DM, E, PDF, N_Events)
    #archiveEnsemblePlot(decay, DM, E, PDF, expectedEvents)

    decay = 'Pi Pi Eta'
    x,y = getPDF(decay, DM)
    f = interp1d(x, y, kind='linear', fill_value='extrapolate')
    E = np.linspace(startE, max(x), 1000)
    PDF = []
    for energy in E:
        PDF.append(f(energy))
    specNorm = integrate.simps(PDF, E)
    expectedEvents = int(round((dRate)/(4*np.pi*DM)*JfactorDracoDecay*specNorm*exposure*secPerYear*0.000404322))
    print('expected events: ' + str(expectedEvents))
    N_Events = np.random.poisson(expectedEvents)
    print(N_Events)

    archiveEnsemblePlot(decay, DM, E, PDF, N_Events)
    #archiveEnsemblePlot(decay, DM, E, PDF, expectedEvents)

    

#Produces plot of simulated photons along with PDF
    '''def gammaPDF(E):
        sig = 1
        return (np.sqrt(2/np.pi)/sig)*np.exp((-E**2/(2*sig**2)))'''
    '''fig, ax = plt.subplots()
    n, bins, patches = ax.hist(E, 100, density=True)
    ax.plot(bins, f(bins)/normalization, '-', label='PDF')
    ax.legend()
    ax.set_title('Simulated Data from PDF')
    ax.set_ylabel('Probability denstity')
    ax.set_xlabel('Energy (Mev)')
    fig.show()'''


#produces plot of etagammaPDF
    '''i=np.linspace(200, 350, 50)
    c=[]
    for val in i:
        c.append(etagammaPDF(val))
    print(integrate.simps(c,i))
    fig, ax = plt.subplots()
    ax.plot(i, c, 'bo')
    ax.set_title('Gamma Distribution w/50 Point Integrand Sample')
    ax.set_ylabel('dN/dE (1/Mev)')
    ax.set_xlabel('E (Mev)')
    ax.legend()
    fig.show()'''
    
#Produces plot of simulated photons along with etagammaPDF
    '''fig, ax = plt.subplots()
    n, bins, patches = ax.hist(randomEnsemble(etagammaPDF, 100), 100, density=True)
    ax.plot(bins, gammaPDF(bins), '--', label='PDF')
    ax.legend()
    ax.set_title('Simulated Data from PDF')
    ax.set_ylabel('Probability denstity')
    ax.set_xlabel('Energy (Mev)')
    fig.show()'''

main()
