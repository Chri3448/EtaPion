import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import csv
from pathlib import Path
import os
import random
random.seed()
import time
from multiprocessing import Pool, Value
import multiprocessing

DM = 820 #<--vestigial
eta = 548
pi = 135
pipm = 140
BRgammagammma = 0.40
BRpipipi = 0.33
BRpi0pipi = 0.23
BRpipi = 1/3


#function f as defined in 2018 paper pg.5 eq(6)
def f(x1, x2, m1, m2, m3, s):
    f_values = np.ones((np.size(x1), np.size(x2)))
    x_values = np.ones((np.size(x1), np.size(x2), 2))
    x_values[:,:,0] = np.tile(x1, (np.size(x2),1)).T
    x_values[:,:,1] = np.tile(x2, (np.size(x1),1))
    f_values = np.where(np.logical_or(s*x_values[:,:,0]**2/4-m1**2 < 0, s*x_values[:,:,1]**2/4-m2**2 < 0),
                        0., f_values)
    heaviside1 = 2-x_values[:,:,0]-x_values[:,:,1]-(2/np.sqrt(s))*np.sqrt(
        m3**2+(np.sqrt(s*x_values[:,:,0]**2/4-m1**2)-np.sqrt(s*x_values[:,:,1]**2/4-m2**2))**2) 
    heaviside2 = x_values[:,:,0]+x_values[:,:,1]+2/np.sqrt(s)*np.sqrt(
        m3**2+(np.sqrt(s*x_values[:,:,0]**2/4-m1**2)+np.sqrt(s*x_values[:,:,1]**2/4-m2**2))**2)-2
    f_values = np.where(np.logical_or(heaviside1 <= 0, heaviside2 <= 0),
                        0., f_values)
    return f_values


#eta injection spectrum for DM->eta pi pi as defined in 2018 paper pg.5 eq(8)
def dNdx1(x1, m1, m2, m3, s, dx1dx2f, M2):
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    i = np.linspace(x2min, x2max, 100)
    c = f(x1, i, m1, m2, m3, s)*M2(x1, i)
    dx2f = integrate.trapz(c, i, axis = 1)
    return dx2f/dx1dx2f

#boosted gamma specrtum from etas in DM->eta pi pi : eta->gamma gamma as defined in 2018 paper pg.5 eq(10)
def M2(x1, x2):
    return 1
#dx1dx2feta is a global variable to save on computaion time
dx1dx2feta1 = 0
dx1dx2feta2 = 0
def etaPDF(x1, DMdecay):
    global dx1dx2feta1
    global dx1dx2feta2
    global DM
    global eta
    global pi
    global pipm
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    if DMdecay == 1:
        m2 = pi #MeV
        m3 = pi #Mev
    elif DMdecay == 2:
        m2 = pipm #MeV
        m3 = pipm #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2feta1 == 0 and DMdecay == 1:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2feta1 = integrate.simps(mArray, list2)
    if dx1dx2feta2 == 0 and DMdecay == 2:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2feta2 = integrate.simps(mArray, list2)
    if DMdecay == 1:
        return dNdx1(x1, m1, m2, m3, s, dx1dx2feta1, M2)
    if DMdecay == 2:
        return dNdx1(x1, m1, m2, m3, s, dx1dx2feta2, M2)
        
def etagammaPDF(E, DMdecay):
    global dx1dx2feta1
    global dx1dx2feta2
    global DM
    global eta
    global pi
    global pipm
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    if DMdecay == 1:
        m2 = pi #MeV
        m3 = pi #Mev
    elif DMdecay == 2:
        m2 = pipm #MeV
        m3 = pipm #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2feta1 == 0 and DMdecay == 1:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2feta1 = integrate.simps(mArray, list2)
    if dx1dx2feta2 == 0 and DMdecay == 2:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2feta2 = integrate.simps(mArray, list2)
    #main integration
    lowerbound = m1/np.sqrt(s)*(2*E/m1+m1/(2*E))
    samplenumber = 1000
    i = np.linspace(lowerbound, x1max, samplenumber)
    if DMdecay == 1:
        c = 2*dNdx1(i, m1, m2, m3, s, dx1dx2feta1, M2)/np.sqrt((s*i**2/4)-m1**2)
        return integrate.trapz(c, i, axis = 0)
    if DMdecay == 2:
        c = 2*dNdx1(i, m1, m2, m3, s, dx1dx2feta2, M2)/np.sqrt((s*i**2/4)-m1**2)
        return integrate.trapz(c, i, axis = 0)

#boosted gamma specrtum from pions in DM->eta pi pi : eta->pi0 pi0 pi0 : pi0 -> gamma gamma
def Mpipipi2(x1, x2):
    global eta
    global pi
    alpha = -0.0318
    x_values = np.ones((np.size(x1), np.size(x2), 2))
    x_values[:,:,0] = np.tile(x1, (np.size(x2),1)).T
    x_values[:,:,1] = np.tile(x2, (np.size(x1),1))
    z = 6/((1-3*pi**2/eta**2)**2)*((x_values[:,:,0]/2-1/3)**2+
                                   (x_values[:,:,1]/2-1/3)**2+
                                   ((2-x_values[:,:,0]-x_values[:,:,1])/2-1/3)**2)
    return 1+2*alpha*z

#dx1dx2fpipipi is a global variable to save on computaion time
dx1dx2fpipipi = 0
#eta rest frame pi distrobution
def RpipipiPDF(x1):
    global dx1dx2fpipipi
    global DM
    global eta
    global pi
    invariantmass = eta #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2fpipipi == 0:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*Mpipipi2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2fpipipi = integrate.simps(mArray, list2)
    #boosting the pion energy spectrum by the eta energy spectrum
    return 3*dNdx1(x1, m1, m2, m3, s, dx1dx2fpipipi, Mpipipi2)

#pipipi spectrum as function of pion energy Ep boosted by a single eta energy x1e
def SBpipipiPDF(Ep, x1e):
    global DM
    global eta
    global pi
    if Ep < pi or Ep == pi:
        return 0
    Pp = np.sqrt(Ep**2-pi**2)
    Ee = DM*x1e/2
    Pe = np.sqrt(Ee**2-eta**2)
    gamma = Ee/eta
    beta = np.sqrt(1-(1/gamma)**2)
    x1max = (2/eta)*gamma*(Ep+Pp*beta)
    x1min = (2/eta)*gamma*(Ep-Pp*beta)
    samplenumber = 50
    i=np.linspace(x1min, x1max, samplenumber)
    c = 0.5*RpipipiPDF(i)/(np.sqrt((eta**2/4)*i**2-pi**2)*gamma*beta)
    return integrate.trapz(c, i, axis = 0)

#pipipi spectrum as function of pion energy Ep boosted by the entire eta energy spectrum
def pipipiPDF(Ep, DMdecay):
    global DM
    global eta
    global pi
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber).tolist()
    #removes zero boost point because of singularities
    i.pop(0)
    i = np.asarray(i)
    c = []
    for val in i:
        c.append(etaPDF(val, DMdecay)*SBpipipiPDF(Ep, val))
    c = np.where(np.isnan(c), 0, c)
    return integrate.simps(c, i, axis = 0)

def pipipigammaPDF(E, DMdecay):
    global DM
    global eta
    global pi
    invariantmass = eta #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    Rx1max = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    REmax = np.sqrt(s)*Rx1max/2
    RPmax = np.sqrt(REmax**2-pi**2)
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    x1maxeta = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    Emaxeta = DM*x1maxeta/2
    gamma = Emaxeta/eta
    beta = np.sqrt(1-(1/gamma)**2)
    upperbound = gamma*(REmax+RPmax*beta)
    lowerbound = (4*E**2+pi**2)/(4*E)
    samplenumber = 15
    i=np.linspace(lowerbound, upperbound, samplenumber).tolist()
    #removes zero boost point because of singularities
    if lowerbound == pi:
        i.pop(0)
    i = np.asarray(i)
    c = []
    for val in i:
        c.append(2*pipipiPDF(val, DMdecay)/np.sqrt(val**2-pi**2))
    return integrate.simps(c, i, axis = 0)[0]

#boosted gamma specrtum from pions in DM->eta pi pi : eta->pi0 pi+ pi- : pi0 -> gamma gamma
def Mpi0pipi2(x1, x2):
    global eta
    global pi
    global pipm
    x_values = np.ones((np.size(x1), np.size(x2), 2))
    x_values[:,:,0] = np.tile(x1, (np.size(x2),1)).T
    x_values[:,:,1] = np.tile(x2, (np.size(x1),1))
    a = -0.94
    b = 0.22
    y = (3/2)*x_values[:,:,0]*eta/(eta-2*pipm-pi)
    return 1+a*y+b*y**2

#dx1dx2fpi0pipi is a global variable to save on computaion time
dx1dx2fpi0pipi = 0
#eta rest frame pi0 distrobution
def Rpi0pipiPDF(x1):
    global dx1dx2fpi0pipi
    global DM
    global eta
    global pi
    global pipm
    invariantmass = eta #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pipm #MeV
    m3 = pipm #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2fpi0pipi == 0:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*Mpi0pipi2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2fpi0pipi = integrate.simps(mArray, list2)
    #boosting the pion energy spectrum by the eta energy spectrum
    return dNdx1(x1, m1, m2, m3, s, dx1dx2fpi0pipi, Mpi0pipi2)

#pi0pipi spectrum as function of pion energy Ep boosted by a single eta energy x1e
def SBpi0pipiPDF(Ep, x1e):
    global DM
    global eta
    global pi
    if Ep < pi or Ep == pi:
        return 0
    Pp = np.sqrt(Ep**2-pi**2)
    Ee = DM*x1e/2
    Pe = np.sqrt(Ee**2-eta**2)
    gamma = Ee/eta
    beta = np.sqrt(1-(1/gamma)**2)
    x1max = (2/eta)*gamma*(Ep+Pp*beta)
    x1min = (2/eta)*gamma*(Ep-Pp*beta)
    samplenumber = 50
    i=np.linspace(x1min, x1max, samplenumber)
    c = 0.5*Rpi0pipiPDF(i)/(np.sqrt((eta**2/4)*i**2-pi**2)*gamma*beta)
    return integrate.trapz(c, i, axis = 0)

#pi0pipi spectrum as function of pion energy Ep boosted by the entire eta energy spectrum
def pi0pipiPDF(Ep, DMdecay):
    global DM
    global eta
    global pi
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber).tolist()
    #removes zero boost point because of singularities
    i.pop(0)
    i = np.asarray(i)
    c = []
    for val in i:
        c.append(etaPDF(val, DMdecay)*SBpi0pipiPDF(Ep, val))
    c = np.where(np.isnan(c), 0, c)
    return integrate.simps(c, i, axis = 0)

def pi0pipigammaPDF(E, DMdecay):
    global DM
    global eta
    global pi
    global pipm
    invariantmass = eta #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pipm #MeV
    m3 = pipm #Mev
    Rx1max = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    REmax = np.sqrt(s)*Rx1max/2
    RPmax = np.sqrt(REmax**2-pi**2)
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = eta #MeV
    m2 = pi #MeV
    m3 = pi #Mev
    x1maxeta = np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2)
    Emaxeta = DM*x1maxeta/2
    gamma = Emaxeta/eta
    beta = np.sqrt(1-(1/gamma)**2)
    upperbound = gamma*(REmax+RPmax*beta)
    lowerbound = (4*E**2+pi**2)/(4*E)
    samplenumber = 15
    i=np.linspace(lowerbound, upperbound, samplenumber).tolist()
    #removes zero boost point because of singularities
    if lowerbound == pi:
        i.pop(0)
    i = np.asarray(i)
    c = []
    for val in i:
        c.append(2*pi0pipiPDF(val, DMdecay)/np.sqrt(val**2-pi**2))
    return integrate.simps(c, i, axis = 0)[0]

#dx1dx2fpipi is a global variable to save on computaion time
dx1dx2fpipi = 0
def pipiPDF(x1):
    global dx1dx2fpipi
    global DM
    global eta
    global pi
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pi #MeV
    m3 = eta #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = 2*(1-(m2+m3)/np.sqrt(s))
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2fpipi == 0:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.simps(c, list1, axis = 0)
        dx1dx2fpipi = integrate.simps(mArray, list2)
    return 2*dNdx1(x1, m1, m2, m3, s, dx1dx2fpipi, M2)

def pipigammaPDF(E):
    global dx1dx2fpipi
    global DM
    global eta
    global pi
    invariantmass = DM #MeV
    s = invariantmass**2
    m1 = pi #MeV
    m2 = pi #MeV
    m3 = eta #Mev
    x1min = 2*m1/np.sqrt(s)
    x1max = 2*(1-(m2+m3)/np.sqrt(s))
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    #preliminary integration to save computation time
    if dx1dx2fpipi == 0:
        list1 = np.linspace(x1min, x1max, 1000)
        list2 = np.linspace(x2min, x2max, 1000)
        c = f(list1, list2, m1, m2, m3, s)*M2(list1, list2)
        mArray = integrate.trapz(c, list1, axis = 0)
        dx1dx2fpipi = integrate.simps(mArray, list2)
    #main integration
    lowerbound = m1/np.sqrt(s)*(2*E/m1+m1/(2*E))
    samplenumber = 1000
    i = np.linspace(lowerbound, x1max, samplenumber)
    c = 2*2*dNdx1(i, m1, m2, m3, s, dx1dx2fpipi, M2)/np.sqrt((s*i**2/4)-m1**2)
    return integrate.trapz(c, i, axis = 0)

#PDF of all eta decays with branching ratios taken in to account, normalized to one photon
def totalPiPiEtagammaPDF(E):
    global BRgammagammma
    global BRpipipi
    global BRpi0pipi
    global BRpipi
    BRtotal = BRgammagammma+BRpipipi+BRpi0pipi
    PDF = BRpipi*(pipigammaPDF(E)+(etagammaPDF(E, 1)*BRgammagammma+pipipigammaPDF(E, 1)*BRpipipi
                 +pi0pipigammaPDF(E, 1)*BRpi0pipi)/BRtotal)+(1-BRpipi)*(etagammaPDF(E, 2)*BRgammagammma
                                                                       +pipipigammaPDF(E, 2)*BRpipipi
                                                                       +pi0pipigammaPDF(E, 2)*BRpi0pipi)/BRtotal
    return PDF

#returns PDF data where each sub-decay is normalized to the proper number of photons. This should improve the
#relative magnitudes of the sub-decays and the overall shape of the total spectrum. Don't start with 0.
def NormalizedTotalGammaPDF(start, stop, step):
    global BRgammagammma
    global BRpipipi
    global BRpi0pipi
    global BRpipi
    BRtotal = BRgammagammma+BRpipipi+BRpi0pipi
    i=np.linspace(start, stop, step)
    ones = []
    twos = []
    for val in i:
        ones.append((val, 1))
        twos.append((val, 2))
    c1=[]
    c2=[]
    c3=[]
    c4=[]
    c5=[]
    c6=[]
    c7=[]
    with Pool(multiprocessing.cpu_count() - 3, resetConstants, (DM,)) as p:
        c1 = p.map(pipigammaPDF, i)
        c2 = p.starmap(etagammaPDF, ones)
        c3 = p.starmap(pipipigammaPDF, ones)
        c4 = p.starmap(pi0pipigammaPDF, ones)
        c5 = p.starmap(etagammaPDF, twos)
        c6 = p.starmap(pipipigammaPDF, twos)
        c7 = p.starmap(pi0pipigammaPDF, twos)
    #print(np.shape(c1),np.shape(c2),np.shape(c3),np.shape(c4),np.shape(c5),np.shape(c6),np.shape(c7))
    #print(c1,c2,c3,c4,c5,c6,c7)
    c1Norm = integrate.simps(c1,i)
    c2Norm = integrate.simps(c2,i)
    c3Norm = integrate.simps(c3,i)
    c4Norm = integrate.simps(c4,i)
    c5Norm = integrate.simps(c5,i)
    c6Norm = integrate.simps(c6,i)
    c7Norm = integrate.simps(c7,i)
    c1 = list(np.array(c1)*(4/c1Norm))
    c2 = list(np.array(c2)*(2/c2Norm))
    c3 = list(np.array(c3)*(6/c3Norm))
    c4 = list(np.array(c4)*(2/c4Norm))
    c5 = list(np.array(c5)*(2/c5Norm))
    c6 = list(np.array(c6)*(6/c6Norm))
    c7 = list(np.array(c7)*(2/c7Norm))
    PDF = list(BRpipi*(np.array(c1)+(np.array(c2)*BRgammagammma+np.array(c3)*BRpipipi+np.array(c4)*BRpi0pipi)/BRtotal)
    +(1-BRpipi)*(np.array(c5)*BRgammagammma+np.array(c6)*BRpipipi+np.array(c7)*BRpi0pipi)/BRtotal)
    return i, PDF
    

def archivePDF(x, PDF):
    global DM
    PDFdir = Path("Pi Pi Eta PDFs/")
    PDFfile = PDFdir / (str(DM) + ".csv")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    with open(PDFfile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(x)
        writer.writerow(PDF)

def archivePlot(x, PDF):
    global DM
    PDFdir = Path("Pi Pi Eta PDFs/")
    PDFfile = PDFdir / (str(DM)+'.png')
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')
    fig, ax = plt.subplots()
    ax.plot(x, PDF, 'k.')
    ax.set_title('DM(' + str(DM) + ' MeV) -> gamma via pi pi eta')
    ax.set_ylabel('dN/dE (1/Mev)')
    ax.set_xlabel('E (Mev)')
    ax.legend()
    fig.savefig(PDFfile)
    #fig.show()

def resetConstants(mass):
    global DM
    global dx1dx2feta1
    global dx1dx2feta2
    global dx1dx2fpipipi
    global dx1dx2fpi0pipi
    global dx1dx2fpipi
    DM = mass
    dx1dx2feta1 = 0
    dx1dx2feta2 = 0
    dx1dx2fpipipi = 0
    dx1dx2fpi0pipi = 0
    dx1dx2fpipi = 0

def combine_pdfs():
    PDFdir = Path("Pi Pi Eta PDFs/")
    E = np.array([])
    Mass = np.array([])
    PDF = np.array([])
    for x in PDFdir.iterdir():
        if str(x).endswith('.csv'):
            mass = float(str(x)[len(str(PDFdir))+1:-4])
            PDFfile = PDFdir / (str(mass) + ".csv")
            e = []
            pdf = []
            with open(PDFfile, "r") as archive:
                reader = csv.reader(archive, delimiter=',')
                for row in reader:
                    if e == []:
                        print(row)
                        for val in row:
                            e.append(float(val))
                        continue
                if pdf == []:
                    for val in row:
                        pdf.append(float(val))
            

    '''
    with open(PDFfile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        reader.readrow()
        reader.readrow()
    '''

def main():
    global DM
    combine_pdfs()
    '''
    DM = 2000
    etaEmax = (DM/2)*np.sqrt(4*eta**2/(DM**2)+(1-(eta**2+2*(pi**2+pi**2))/(DM**2))**2-8*eta**2*(pi**2+pi**2)/(DM**2)**2)
    etaBeta = np.sqrt(((etaEmax/eta)**(-2)-1)*(-1))
    photonEmax = (eta/2)*np.sqrt((1+etaBeta)/(1-etaBeta))
    print(time.asctime(time.localtime()))
    i = np.linspace(0.01, photonEmax + 20, 50)
    c = []
    for val in i:
        c.append(totalPiPiEtagammaPDF(val))
    BRtotal = BRgammagammma+BRpipipi+BRpi0pipi
    print('average number of photons: ' + str((2*BRgammagammma+6*BRpipipi+2*BRpi0pipi)/BRtotal + BRpipi*4))
    print('integral of PDF: ' + str(integrate.simps(c,i)))
    print(time.asctime(time.localtime()))
    '''
    '''
    if __name__=='__main__':
        print(time.asctime(time.localtime()))
        for mass in np.linspace(880, 880, 1):
            DM = mass
            resetConstants(mass)
            print(DM)
            etaEmax = (DM/2)*np.sqrt(4*eta**2/(DM**2)+(1-(eta**2+2*(pi**2+pi**2))/(DM**2))**2-8*eta**2*(pi**2+pi**2)/(DM**2)**2)
            etaBeta = np.sqrt(((etaEmax/eta)**(-2)-1)*(-1))
            photonEmax = (eta/2)*np.sqrt((1+etaBeta)/(1-etaBeta))
            i, c = NormalizedTotalGammaPDF(0.01, photonEmax + 20, 50)
            BRtotal = BRgammagammma+BRpipipi+BRpi0pipi
            print('average number of photons: ' + str((2*BRgammagammma+6*BRpipipi+2*BRpi0pipi)/BRtotal + BRpipi*4))
            print('integral of PDF: ' + str(integrate.simps(c,i)))
            archivePDF(i,c)
            archivePlot(i,c)
        print(time.asctime(time.localtime()))
    '''

#Produces plot of simulated photons along with PDF
    '''def gammaPDF(E):
        sig = 1
        return (np.sqrt(2/np.pi)/sig)*np.exp((-E**2/(2*sig**2)))

    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(randomEnsemble(gammaPDF, 100000), 100, density=True)
    ax.plot(bins, gammaPDF(bins), '--', label='PDF')
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
