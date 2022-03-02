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
    if s*x1**2/4-m1**2 < 0 or s*x2**2/4-m2**2 < 0:
        return 0
    heaviside1 = 2-x1-x2-(2/np.sqrt(s))*np.sqrt(m3**2+(np.sqrt(s*x1**2/4-m1**2)-np.sqrt(s*x2**2/4-m2**2))**2) 
    heaviside2 = x1+x2+2/np.sqrt(s)*np.sqrt(m3**2+(np.sqrt(s*x1**2/4-m1**2)+np.sqrt(s*x2**2/4-m2**2))**2)-2
    if heaviside1 > 0 and heaviside2 > 0:
        return 1
    else:
        return 0
#Produces plot of f
'''print(f(2*(549)/np.sqrt(s), 2*(141)/np.sqrt(s), m1, m2, m3, s))
i=np.linspace(2*m1/np.sqrt(s), np.sqrt(4*m1**2/s+(1-(m1**2+4*m2**2)/s)**2-16*m1**2*m2**2/s**2), 1000)
c=[]
for val in i:
    c.append(f(val, 2*(144)/np.sqrt(s), m1, m2, m3, s))
fig, ax = plt.subplots()
ax.plot(i, c, 'ro')
ax.legend()
fig.show()'''

#eta injection spectrum for DM->eta pi pi as defined in 2018 paper pg.5 eq(8)
def dNdx1(x1, m1, m2, m3, s, dx1dx2f, M2):
    x1min = 2*m1/np.sqrt(s)
    x1max = np.sqrt(4*m1**2/s+(1-(m1**2+2*(m2**2+m3**2))/s)**2-8*m1**2*(m2**2+m3**2)/s**2)
    x2min = 2*m2/np.sqrt(s)
    x2max = 2*(1-(m1+m3)/np.sqrt(s))
    '''dx2f = integrate.quad(lambda x2: f(x1, x2, m1, m2, m3, s)*M2(x1, x2), x2min, x2max)[0]'''
    samplenumber = 100
    i=np.linspace(x2min, x2max, samplenumber)
    c=[]
    for val in i:
        c.append(f(x1, val, m1, m2, m3, s)*M2(x1, val))
    dx2f = integrate.simps(c, i)
    return dx2f/dx1dx2f
#finds endpoints such that dNdx1 is smooth between them
'''
endpt1 = 0
endpt2 = 2
i = np.linspace(x1min, x1min + .001, 1000)
for val in i:
    if not dNdx1(val, m1, m2, m3, s, dx1dx2f, M2) == 0:
        endpt1 = val
        break
i = np.linspace(x1max - .001, x1max, 1000)
for val in i:
    if dNdx1(val, m1, m2, m3, s, dx1dx2f, M2) == 0:
        break
    endpt2 = val
print(endpt1, x1min, endpt2, x1max)'''

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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2feta1 = integrate.simps(mArray, list2)
    if dx1dx2feta2 == 0 and DMdecay == 2:
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2feta2 = integrate.simps(mArray, list2)
        '''dx1dx2feta = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*M2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2feta calculated')
        print(dx1dx2feta)'''
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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2feta1 = integrate.simps(mArray, list2)
    if dx1dx2feta2 == 0 and DMdecay == 2:
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2feta2 = integrate.simps(mArray, list2)
        '''dx1dx2feta = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*M2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2feta calculated')
        print(dx1dx2feta)'''
    #main integration
    lowerbound = m1/np.sqrt(s)*(2*E/m1+m1/(2*E))
    samplenumber = 100
    i=np.linspace(lowerbound, x1max, samplenumber)
    c=[]
    if DMdecay == 1:
        for val in i:
            c.append(2*dNdx1(val, m1, m2, m3, s, dx1dx2feta1, M2)/np.sqrt((s*val**2/4)-m1**2))
        return integrate.simps(c, i)
    if DMdecay == 2:
        for val in i:
            c.append(2*dNdx1(val, m1, m2, m3, s, dx1dx2feta2, M2)/np.sqrt((s*val**2/4)-m1**2))
        return integrate.simps(c, i)
    '''return 2*integrate.quad(lambda x1: dNdx1(x1)/np.sqrt((s*x1**2/4)-m1**2),
                        lowerbound, endpt)'''

#boosted gamma specrtum from pions in DM->eta pi pi : eta->pi0 pi0 pi0 : pi0 -> gamma gamma
def Mpipipi2(x1, x2):
    global eta
    global pi
    alpha = -0.0318
    z = 6/((1-3*pi**2/eta**2)**2)*((x1/2-1/3)**2+(x2/2-1/3)**2+((2-x1-x2)/2-1/3)**2)
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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*Mpipipi2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2fpipipi = integrate.simps(mArray, list2)
        '''dx1dx2fpipipi = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*Mpipipi2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2fpipipi calculated')
        print(dx1dx2fpipipi)'''
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
    '''return integrate.quad(lambda x1: 0.5*RpipipiPDF(x1)/(np.sqrt((eta**2/4)*x1**2-pi**2)*gamma*beta), x1min, x1max)[0]'''
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber)
    c=[]
    for val in i:
        PDFval = RpipipiPDF(val)
        if PDFval == 0:
            c.append(0)
        else:
            c.append(0.5*PDFval/(np.sqrt((eta**2/4)*val**2-pi**2)*gamma*beta))
    return integrate.simps(c, i)

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
    '''return integrate.quad(lambda x1: etaPDF(x1)*SBpipipiPDF(Ep, x1), x1min, x1max)[0]'''
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber).tolist()
    #removes zero boost point because of singularities
    i.pop(0)
    np.asarray(i)
    c=[]
    for val in i:
        PDFval1 = etaPDF(val, DMdecay)
        PDFval2 = SBpipipiPDF(Ep, val)
        if PDFval1 == 0 or np.isnan(PDFval1) or PDFval2 == 0 or np.isnan(PDFval2):
            c.append(0)
        else:
            c.append(PDFval1*PDFval2)
    return integrate.simps(c, i)

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
    np.asarray(i)
    c=[]
    for val in i:
        c.append(2*pipipiPDF(val, DMdecay)/np.sqrt(val**2-pi**2))
    return integrate.simps(c, i)

#boosted gamma specrtum from pions in DM->eta pi pi : eta->pi0 pi+ pi- : pi0 -> gamma gamma
def Mpi0pipi2(x1, x2):
    global eta
    global pi
    global pipm
    a = -0.94
    b = 0.22
    y = (3/2)*x1*eta/(eta-2*pipm-pi)
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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*Mpi0pipi2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2fpi0pipi = integrate.simps(mArray, list2)
        '''dx1dx2fpi0pipi = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*Mpi0pipi2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2fpi0pipi calculated')
        print(dx1dx2fpi0pipi)'''
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
    '''return integrate.quad(lambda x1: 0.5*Rpi0pipiPDF(x1)/(np.sqrt((eta**2/4)*x1**2-pi**2)*gamma*beta), x1min, x1max)[0]'''
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber)
    c=[]
    for val in i:
        PDFval = Rpi0pipiPDF(val)
        if PDFval == 0:
            c.append(0)
        else:
            c.append(0.5*PDFval/(np.sqrt((eta**2/4)*val**2-pi**2)*gamma*beta))
    return integrate.simps(c, i)

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
    '''return integrate.quad(lambda x1: etaPDF(x1)*SBpi0pipiPDF(Ep, x1), x1min, x1max)[0]'''
    samplenumber = 15
    i=np.linspace(x1min, x1max, samplenumber).tolist()
    #removes zero boost point because of singularities
    i.pop(0)
    np.asarray(i)
    c=[]
    for val in i:
        PDFval1 = etaPDF(val, DMdecay)
        PDFval2 = SBpi0pipiPDF(Ep, val)
        if PDFval1 == 0 or np.isnan(PDFval1) or PDFval2 == 0 or np.isnan(PDFval2):
            c.append(0)
        else:
            c.append(PDFval1*PDFval2)
    return integrate.simps(c, i)

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
    np.asarray(i)
    c=[]
    for val in i:
        c.append(2*pi0pipiPDF(val, DMdecay)/np.sqrt(val**2-pi**2))
    return integrate.simps(c, i)

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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2fpipi = integrate.simps(mArray, list2)
        '''dx1dx2fpipi = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*M2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2fpipi calculated')
        print(dx1dx2fpipi)'''
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
        list1 = np.linspace(x1min, x1max, 500)
        list2 = np.linspace(x2min, x2max, 500)
        mArray = []
        for a in list2:
            c = []
            for b in list1:
                c.append(f(b, a, m1, m2, m3, s)*M2(b, a))
            mArray.append(integrate.simps(c, list1))
        dx1dx2fpipi = integrate.simps(mArray, list2)
        '''dx1dx2fpipi = integrate.dblquad(lambda x_1, x2: f(x_1, x2, m1, m2, m3, s)*M2(x_1, x2), x2min, x2max,
                                       lambda x2: x1min,lambda x2: x1max)[0]
        print('dx1dx2fpipi calculated')
        print(dx1dx2fpipi)'''
    #main integration
    lowerbound = m1/np.sqrt(s)*(2*E/m1+m1/(2*E))
    samplenumber = 15
    i=np.linspace(lowerbound, x1max, samplenumber)
    c=[]
    for val in i:
        c.append(2*2*dNdx1(val, m1, m2, m3, s, dx1dx2fpipi, M2)/np.sqrt((s*val**2/4)-m1**2))
    return integrate.simps(c, i)
    '''return 2*integrate.quad(lambda x1: dNdx1(x1)/np.sqrt((s*x1**2/4)-m1**2),
                        lowerbound, endpt)'''

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
    with Pool(9, resetConstants, (DM,)) as p:
        c1 = p.map(pipigammaPDF, i)
        c2 = p.starmap(etagammaPDF, ones)
        c3 = p.starmap(pipipigammaPDF, ones)
        c4 = p.starmap(pi0pipigammaPDF, ones)
        c5 = p.starmap(etagammaPDF, twos)
        c6 = p.starmap(pipipigammaPDF, twos)
        c7 = p.starmap(pi0pipigammaPDF, twos)
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
    PDFfile = PDFdir / (str(DM) + ".png")
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

def main():
    global DM
    if __name__=='__main__':
        print(time.asctime(time.localtime()))
        for mass in np.linspace(830, 1100, 28):
            DM = mass
            resetConstants(mass)
            print(DM)
            etaEmax = (DM/2)*np.sqrt(4*eta**2/(DM**2)+(1-(eta**2+2*(pi**2+pi**2))/(DM**2))**2-8*eta**2*(pi**2+pi**2)/(DM**2)**2)
            etaBeta = np.sqrt(((etaEmax/eta)**(-2)-1)*(-1))
            photonEmax = (eta/2)*np.sqrt((1+etaBeta)/(1-etaBeta))
            #i, c = NormalizedTotalGammaPDF(0.01, photonEmax + 20, 50)
            i, c = NormalizedTotalGammaPDF(0.01, 1000, 500)
            BRtotal = BRgammagammma+BRpipipi+BRpi0pipi
            print('average number of photons: ' + str((2*BRgammagammma+6*BRpipipi+2*BRpi0pipi)/BRtotal + BRpipi*4))
            print('integral of PDF: ' + str(integrate.simps(c,i)))
            archivePDF(i,c)
            archivePlot(i,c)
        print(time.asctime(time.localtime()))
    

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
