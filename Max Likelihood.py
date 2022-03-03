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
    return x,y

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

def smearBackground(energy, E, PDF, res):
    sig = res/(2*np.sqrt(2*np.log(2)))
    f = interp1d(E, PDF, kind='linear', fill_value='extrapolate')
    E = np.linspace(startE, max(E), 1000)
    PDF = f(E)/np.sqrt(2*np.pi*(sig*E)**2)*np.exp(-.5*(energy - E)**2/(sig*E)**2)
    return integrate.trapz(PDF, E)

def getSmearedBackground(exposure, Sangle, res):
    sig = res/(2*np.sqrt(2*np.log(2)))
    E = np.linspace(startE/(1+4*sig), 1000, 1000)
    PDF = background(E, exposure, Sangle)
    E2 = np.linspace(startE, max(E) + 4*sig*max(E), 1000)
    PDF2 = smearPDF(E2, E, PDF, res)
    return E2, PDF2

def smearEnsemble(ensemble, res):
    sig = res/(2*np.sqrt(2*np.log(2)))
    newEnsemble = np.random.normal(ensemble, sig*ensemble)
    indices = np.where(newEnsemble >= startE)
    newEnsemble = newEnsemble[indices]
    #print('Photons lost in smear: ' + str(len(ensemble) - len(newEnsemble)))
    return newEnsemble

def getMasses(decay):
    PDFdir = Path(decay + " PDFs/")
    masses = []
    for x in PDFdir.iterdir():
        if str(x).endswith('.csv'):
            mass = float(str(x)[len(str(PDFdir))+1:-4])
            if not mass > 1100.0:
                masses.append(mass)
    return np.sort(np.array(masses))

def draw_from_pdf(Pc, Ndraws):
        # draw random counts from P(c)
        cdf = np.cumsum(Pc)
        rands = np.random.rand(Ndraws)
        # Draw Ndraws times from Pc
        d_vec = np.searchsorted(cdf, rands)
        return d_vec

def getEnsemble(decay, DM, dRate, exposure, Sangle):
    '''
    Ensemble = []
    Edir = Path(decay + " Ensembles/")
    Efile = Edir / (str(DM) + ".csv")
    with open(Efile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        for row in reader:
            if Ensemble == []:
                for val in row:
                    Ensemble.append(float(val))
    '''
    '''
    secPerYear = 3.154*10**7
    JfactorDracoDecay = 5.77*(10**24) #(MeV cm^-2 sr^-1)
    x,y = getPDF(decay, DM)
    f = interp1d(x, y, kind='linear', fill_value='extrapolate')
    E = np.linspace(startE, max(x), 1000)
    PDF = f(E)
    specNorm = integrate.simps(PDF, E)
    expectedEvents = int(round((dRate)/(4*np.pi*DM)*JfactorDracoDecay*specNorm*exposure*secPerYear*Sangle))
    N_Events = np.random.poisson(expectedEvents)
    Edir = Path()
    Efile = Edir / (decay+' Ensembles') / (str(DM)+".csv")
    with open(Efile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        for row in reader:
            raw = np.asarray(row)
            break
    indices = np.random.randint(low = 0, high = len(raw) - 1, size = N_Events)
    Ensemble = raw[indices].astype('float')
    return Ensemble
    '''
    secPerYear = 3.154*10**7
    JfactorDracoDecay = 5.77*(10**24) #(MeV cm^-2 sr^-1)
    x,y = getPDF(decay, DM)
    f = interp1d(x, y, kind='linear', fill_value='extrapolate')
    E = np.linspace(min(x), max(x), 1000)
    PDF = f(E)
    specNorm = integrate.simps(PDF, E)
    expectedEvents = int(round((dRate)/(4*np.pi*DM)*JfactorDracoDecay*specNorm*exposure*secPerYear*Sangle))
    N_Events = np.random.poisson(expectedEvents)
    Edir = Path(decay + " Ensembles/")
    Efile = Edir / (str(DM) + ".csv")
    Eplotfile = Edir / (str(DM) + ".png")
    '''
    try:
        os.mkdir(Edir)
    except:
        print(Edir,'already exists')
    '''
    newE = np.linspace(min(E), max(E), 100000)
    newPDF = f(newE)
    Ensemble = newE[draw_from_pdf(newPDF/np.sum(newPDF), N_Events)]
    '''
    with open(Efile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(Ensemble)
    '''
    '''
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(Ensemble, 1000, density=True, label='Photon Histogram')
    normalization = integrate.simps(newPDF, newE)
    ax.plot(bins, f(bins)/normalization, '-', label='PDF')
    ax.legend()
    ax.set_title(str(N_Events) + ' Simulated Photons from DM(' + str(DM) + ' MeV) -> gamma via ' + decay)
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlabel('Energy (Mev)')
    fig.savefig(Eplotfile)
    fig.show()
    '''
    return Ensemble


def getBackEnsemble(exposure, Sangle):
    '''
    Ensemble = []
    Edir = Path()
    Efile = Edir / ("smeared background.csv")
    with open(Efile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        for row in reader:
            if Ensemble == []:
                for val in row:
                    Ensemble.append(float(val))
    '''
    E = np.linspace(startE/2, 1000, 10000)
    PDF = background(E, exposure, Sangle)
    specNorm = integrate.simps(PDF, E)
    #print("Background Integral: " + str(specNorm))
    expectedEvents = specNorm
    N_Events = np.random.poisson(expectedEvents)
    Ensemble = E[draw_from_pdf(PDF/np.sum(PDF), N_Events)]
    Edir = Path()
    Efile = Edir / ("background.csv")
    Eplotfile = Edir / ("background.png")
    '''
    with open(Efile, "w+") as archive:
        writer = csv.writer(archive, delimiter=',')
        writer.writerow(Ensemble)
    '''
    '''
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(Ensemble, 1000, density=True, label='Photon Histogram')
    ax.plot(bins, background(bins, exposure, Sangle)/specNorm, '-', label='PDF')
    ax.legend()
    ax.set_title(str(N_Events) + ' Simulated Photons from background')
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlabel('Energy (Mev)')
    fig.show()
    '''
    '''
    with open(Efile, "r") as archive:
        reader = csv.reader(archive, delimiter=',')
        for row in reader:
            raw = np.asarray(row)
            break
    indices = np.random.randint(low = 0, high = len(raw) - 1, size = N_Events)
    Ensemble = raw[indices].astype('float')
    '''
    return Ensemble

#returns the natural log of joint prob dist for just a sungle decay and no background
def JPDsingle(decay, DM, back, Ensemble):         
    E, PDF = getPDF(decay, DM)
    f = interp1d(E, PDF, kind='linear')
    normalization = integrate.simps(PDF, E)
    val = 0
    for photon in Ensemble:
        try:
            prob = f(photon)
            if prob < 0:
                prob = 0
        except:
            val = 0
            break
        if prob == 0:
            val = 0
            break
        val += np.log(prob/normalization)
    return val

#returns the natural log of joint prob dist of model + background
def JPD_2decays(E1, PDF1, E2, PDF2, E3, PDF3, decayRate1, decayRate2, DM, exposure, Sangle, Ensemble):
    JfactorDracoDecay = 5.77*(10**24) #(MeV cm^-2 sr^-1)
    secPerYear = 3.154*10**7
#get spectra:
    '''
    E1, PDF1 = getPDF('Pi Eta', DM)
    E2, PDF2 = getPDF('Pi Pi Eta', DM)
    E3 = np.linspace(startE, 1000, 10000)
    PDF3 = background(E3, exposure, Sangle)
    
    E1, PDF1 = getPDF('Gauss3', DM)
    E2, PDF2 = getPDF('Gauss4', DM)
    '''
    '''
    E1, PDF1 = getSmearedPDF(decay_channels[0], DM, .3)
    E2, PDF2 = getSmearedPDF(decay_channels[1], DM, .3)
    E3, PDF3 = getSmearedBackground(exposure, Sangle, .3)
    '''
    
#get interpolation functons of decay spectra
    f1 = interp1d(E1, PDF1, kind='linear', fill_value='extrapolate')
    f2 = interp1d(E2, PDF2, kind='linear', fill_value='extrapolate')
#get spectral norms:
    '''
    E1 = np.linspace(startE, max(E1), 1000)
    E2 = np.linspace(startE, max(E2), 1000)
    PDF1 = f1(E1)
    PDF2 = f2(E2)
    '''
    specNorm1 = integrate.simps(PDF1, E1)
    specNorm2 = integrate.simps(PDF2, E2)
    specNorm3 = integrate.simps(PDF3, E3)
#get expected number of events:
    expectedEvents1 = decayRate1/(4*np.pi*DM)*JfactorDracoDecay*specNorm1*exposure*secPerYear*Sangle
    expectedEvents2 = decayRate2/(4*np.pi*DM)*JfactorDracoDecay*specNorm2*exposure*secPerYear*Sangle
    expectedEvents3 = specNorm3
    totalExpectedEvents = expectedEvents1 + expectedEvents2 + expectedEvents3
    #totalExpectedEvents = expectedEvents1 + expectedEvents2
    #print(totalExpectedEvents)
#combines background and decays into one PDF based on relative numbers of events
    Emax = max(max(E1), max(E2), max(E3))
    Emin = min(min(E1), min(E2), min(E3))
    if Emin < startE:
        Emin = startE

    E = np.linspace(Emin, Emax, 5000)
    '''
    PDF = []
    for energy in E:
        funcValue1 = 0 #necessary because the extrapolated interpolations and background
        funcValue2 = 0 #are unreliable outside domain
        funcValue3 = 0
        if energy >= min(E1) and energy <= max(E1):
            funcValue1 = f1(energy)
        if energy >= min(E2) and energy <= max(E2):
            funcValue2 = f2(energy)
        if energy >= min(E3) and energy <= max(E3):
            funcValue3 = background(energy, exposure, Sangle)
        PDF.append(((expectedEvents1*funcValue1/specNorm1)
                    +(expectedEvents2*funcValue2/specNorm2)
                    +(expectedEvents3*funcValue3/specNorm3))/totalExpectedEvents)
        #            )/totalExpectedEvents)
    '''
    PDF = ((expectedEvents1*np.where(E > max(E1), 0, f1(E))/specNorm1)
           +(expectedEvents2*np.where(E > max(E2), 0, f2(E))/specNorm2)
           +(expectedEvents3*background(E, exposure, Sangle)/specNorm3))/totalExpectedEvents
#the simpsons integrations give close but not perfect PDFs, and normalization is necessary
    normalization = integrate.simps(PDF, E)
    f = interp1d(E, PDF, kind='linear', fill_value= 'extrapolate')
    #print(expectedEvents1, expectedEvents2, expectedEvents3)
    #print(normalization)
    '''
    val = 0
    for photon in Ensemble:
        try:
            prob = f(photon)
            if prob < 0:
                print('Error 1')
                prob = 0
        except:
            print('Error 2')
            val = 0
            break
        if prob == 0:
            print('Error 3')
            #val = 0
            #break
            continue
        val += np.log(prob/normalization)
    '''
    probs = f(Ensemble)/normalization
    if any(probs <= 0):
        print('Warning: some photons have zoro or negative probability')
    val = np.sum(np.log(probs))
#the true probability is the product of the spectral probabilities and the poisson probability
    #'''
    scipy.stats.poisson
    numberEvents = len(Ensemble)
    logPoisonFactor1 = numberEvents*np.log(totalExpectedEvents)
    logPoisonFactor2 = -totalExpectedEvents
    '''
    logPoisonFactor3 = 0
    for i in range(numberEvents):
        logPoisonFactor3 += np.log(i+1)
    '''
    logPoisonFactor3 = np.sum(np.log(np.arange(numberEvents)+1))
    logPoissonProb = logPoisonFactor1 + logPoisonFactor2 - logPoisonFactor3
    val += logPoissonProb
    #'''

    return val

def getArchive(decay):
    PDFdir = Path(decay + " PDFs/")
    files = os.listdir(PDFdir)
    for file in files:
        if not file.endswith('.csv'):
            files.pop(files.index(file))
    for file in files:
        files[files.index(file)] = float(files[files.index(file)][:-4])
    files.sort()
    return files

#takes function values and inputs as lists, returns maximum of interpolation function and input
def maximum(x, val):
    valArray = -np.array(val)
    func = interp1d(x, valArray, kind='linear',
                      fill_value='extrapolate')
    testPoints = np.linspace(x[5], x[len(x)-1], 10)
    for point in testPoints:
        res = minimize(func, point, method='nelder-mead',
               options={'disp': True})
        print(res)
        if res.success == True:
            break
    if res.success == False:
        return False
    return res.x[0]

#background PDF as a function of energy
#energy (MeV), exposure (cm^2 yr), Sangle (sr) 
def background(energy, exposure, Sangle):
    secPerYear = 3.154*10**7
    preFactor = 0.00274
    back = preFactor*exposure*Sangle*secPerYear/energy**2
    return (1)*back

#gives the probability of getting something at least as extreme in %
#takes log likelihood of most likely point in parameter space and true point
def sigmaProb(mostLikely, true, deg = 2.):
    lam = -2*(true - mostLikely)
    chi2CDF = scipy.special.gammainc(deg/2, lam/2)
    return 1 - chi2CDF

###################################################################################################################
#various plot functions below here
###################################################################################################################

def archiveMaxPlot(DMassPi, PiLikelihood, PiFunc,
                   DMassPiPi, PiPiLikelihood, PiPiFunc):
    '''PDFdir = Path("Max Likelihood Plots/")
    PDFfile = PDFdir / (str(DM) + ".png")
    try:
        os.mkdir(PDFdir)
    except:
        print(PDFdir,'already exists')'''
    fig, ax = plt.subplots()
    ax.plot(DMassPi, PiLikelihood, 'k.', label='DM to Pi Eta')
    ax.plot(DMassPi, PiFunc(DMassPi), '-.', label='DM to Pi Eta')
    ax.plot(DMassPiPi, PiPiLikelihood, 'b.', label='DM to Pi Pi Eta')
    ax.plot(DMassPiPi, PiPiFunc(DMassPiPi), '-.', label='DM to Pi Pi Eta')
    ax.set_title('Log Likelyhood vs Dark Matter Mass')
    ax.set_ylabel('Log Likelihood')
    ax.set_xlabel('E (Mev)')
    ax.legend()
    #fig.savefig(PDFfile)
    fig.show()

def likelihoodPlot_2decays(true_decayRate1, true_decayRate2, decayRates1, decayRates2, DM, exposure, Sangle, resolution):
    Ensemble1 = getEnsemble('Pi Eta', DM, true_decayRate1, exposure, Sangle)
    Ensemble2 = getEnsemble('Pi Pi Eta', DM, true_decayRate2, exposure, Sangle)
    Ensemble3 = getBackEnsemble(exposure, Sangle)
    Ensemble = np.concatenate((Ensemble1, Ensemble2, Ensemble3))
    Ensemble = smearEnsemble(Ensemble, resolution)
    print(str(len(Ensemble1))+', '+str(len(Ensemble2))+', '+str(len(Ensemble3)))
    likelihoods = []
    for decay1 in decayRates1:
        for decay2 in decayRates2:
            print('model coordinates: (' + str(decay1) + ',' + str(decay2) + ')')
            Jprob = JPD_2decays(decay1, decay2, DM, exposure, 0.000404322, Ensemble)
            likelihoods.append([decay1, decay2, Jprob])
    maxJprob = -np.inf
    bestModel = []
    for model in likelihoods:
        print(model)
        if model[2] > maxJprob:
            maxJprob = model[2]
            bestModel = model
    trueModel = [true_decayRate1, true_decayRate2, JPD_2decays(true_decayRate1, true_decayRate2, DM, exposure, Sangle, Ensemble)]
    print('Best model:')
    print(bestModel)
    print('True model:')
    print(trueModel)
    sigma = 1+scipy.special.erf(-(5)/np.sqrt(2))#5 sigma ellipse
    maxJprob = bestModel[2]
    fig, ax = plt.subplots()
    firstBlue = True
    firstBlack = True
    for model in likelihoods:
        if model[2] == bestModel[2]:
            ax.plot(model[0], model[1], 'r.', label='Max Probability')
        elif sigmaProb(bestModel[2], model[2]) > sigma:
            if firstBlue:
                ax.plot(model[0], model[1], 'b.', label='Probability > ' + str(sigma))
                firstBlue = False
            else:
                ax.plot(model[0], model[1], 'b.')
        elif sigmaProb(bestModel[2], model[2]) <= sigma:
            if firstBlack:
                ax.plot(model[0], model[1], 'k.', label='Probability <= ' + str(sigma))
                firstBlack = False
            else:
                ax.plot(model[0], model[1], 'k.')
    ax.set_title('Probability of True Model: {:.5%}'.format(sigmaProb(bestModel[2], trueModel[2]))+'\n'
                 + 'True Model: ({:.5} , {:.5})'.format(float(trueModel[0]), float(trueModel[1]))+'\n'
                 + 'Best Model: ({:.5} , {:.5})'.format(float(bestModel[0]), float(bestModel[1])))
    ax.set_ylabel('Pi Pi Eta Decay Rate [s^-1]')
    ax.set_xlabel('Pi Eta Decay Rate [s^-1]')
    ax.legend(loc='upper right')
    fig.show()
    modelPlot(trueModel, bestModel, DM, exposure, 0.000404322, Ensemble)

def likelihoodPlot_2decays_DMs(decay_channels, true_params, radius, num_points, exposure, Sangle, resolution, marginalize = True):
    DM = true_params[0]
    true_decayRate1 = true_params[1]
    true_decayRate2 = true_params[2]
    decayRates1 = np.linspace(0., true_decayRate1 + 2*radius, num_points)
    decayRates2 = np.linspace(0., true_decayRate2 + 2*radius, num_points)
    DM_thresholds = np.zeros(np.size(decay_channels))
    for di in np.arange(np.size(decay_channels)):
        DM_thresholds[di] = np.min(getMasses(decay_channels[di]))
    DMs = getMasses(decay_channels[np.where(DM_thresholds == np.max(DM_thresholds))[0][0]])
    Ensemble1 = getEnsemble(decay_channels[0], DM, true_decayRate1, exposure, Sangle)
    Ensemble2 = getEnsemble(decay_channels[1], DM, true_decayRate2, exposure, Sangle)
    Ensemble3 = getBackEnsemble(exposure, Sangle)
    Ensemble = np.concatenate((Ensemble1, Ensemble2, Ensemble3))
    Ensemble = smearEnsemble(Ensemble, resolution)
    print('Ensembles: '+ str(len(Ensemble1))+', '+str(len(Ensemble2))+', '+str(len(Ensemble3)))
    likelihoods = []
    for decay1 in decayRates1:
        for decay2 in decayRates2:
            vec = np.array([decay1-true_decayRate1, decay2-true_decayRate2])
            if np.absolute(np.sum(vec*np.array([1,-1]))/np.sqrt(2)) < 2*radius and np.absolute(np.sum(vec*np.array([1,1]))/np.sqrt(2)) < radius/2:
            #if (decay1-true_decayRate1)**2 + (decay2-true_decayRate2)**2 < radius**2:
                mass_likelihoods = []
                for mass in DMs:
                    if True:
                        print('model coordinates: (' + str(decay1) + ',' + str(decay2) + ',' + str(mass) + ')')
                        Jprob = JPD_2decays(decay1, decay2, mass, exposure, Sangle, Ensemble)
                        if marginalize:
                            mass_likelihoods.append(Jprob)
                        else:
                            likelihoods.append([decay1, decay2, mass, Jprob])
                if marginalize:
                    marginalized_Jprob = np.max(mass_likelihoods)+np.log(np.sum(np.exp(np.array(mass_likelihoods)-np.max(mass_likelihoods))))
                    likelihoods.append([decay1, decay2, np.nan, marginalized_Jprob])
    likelihoods = np.array(likelihoods)       
    maxJprob = -np.inf
    bestModel = []
    for model in likelihoods:
        print(model)
        if model[3] > maxJprob:
            maxJprob = model[3]
            bestModel = model
    if marginalize:
        mass_likelihoods = []
        for mass in DMs:
            print('calculating true model for mass ' + str(mass) + 'MeV')
            if True:
                Jprob = JPD_2decays(true_decayRate1, true_decayRate2, mass, exposure, Sangle, Ensemble)
                mass_likelihoods.append(Jprob)
        marginalized_Jprob = np.max(mass_likelihoods)+np.log(np.sum(np.exp(np.array(mass_likelihoods)-np.max(mass_likelihoods))))
        trueModel = [true_decayRate1, true_decayRate2, DM, marginalized_Jprob]
    else:
        trueModel = [true_decayRate1, true_decayRate2, DM, JPD_2decays(true_decayRate1, true_decayRate2, DM, exposure, Sangle, Ensemble)]
    print('Best model:')
    print(bestModel)
    print('True model:')
    print(trueModel)
    sigma = 1+scipy.special.erf(-(5)/np.sqrt(2))#5 sigma ellipse
    maxJprob = bestModel[3]
    fig, ax = plt.subplots()
    firstBlue = True
    firstBlack = True
    for model in likelihoods:
        if model[3] == bestModel[3]:
            ax.plot(model[0], model[1], 'r.', label='Max Probability')
        elif sigmaProb(bestModel[3], model[3]) > sigma:
            if firstBlue:
                ax.plot(model[0], model[1], 'b.', label='Probability > ' + str(sigma))
                firstBlue = False
            else:
                ax.plot(model[0], model[1], 'b.')
        elif sigmaProb(bestModel[3], model[3]) <= sigma:
            if firstBlack:
                ax.plot(model[0], model[1], 'k.', label='Probability <= ' + str(sigma))
                firstBlack = False
            else:
                ax.plot(model[0], model[1], 'k.')
    for decay1 in decayRates1:
        for decay2 in decayRates2:
            vec = np.array([decay1-true_decayRate1, decay2-true_decayRate2])
            if not np.absolute(np.sum(vec*np.array([1,-1]))/np.sqrt(2)) < 2*radius or not np.absolute(np.sum(vec*np.array([1,1]))/np.sqrt(2)) < radius/2:
            #if not (decay1-true_decayRate1)**2 + (decay2-true_decayRate2)**2 < radius**2:
                ax.plot(decay1, decay2, 'g.')
    ax.set_title('Probability of True Model: {:.5%}'.format(sigmaProb(bestModel[3], trueModel[3]))+'\n'
                 + 'True Model: ({:.5} , {:.5}) {:.5}'.format(float(trueModel[0]), float(trueModel[1]), float(DM))+'MeV\n'
                 + 'Best Model: ({:.5} , {:.5})'.format(float(bestModel[0]), float(bestModel[1])))
    ax.set_ylabel('Pi Pi Eta Decay Rate [s^-1]')
    ax.set_xlabel('Pi Eta Decay Rate [s^-1]')
    #ax.legend(loc='lower left')
    fig.show()
    #modelPlot(trueModel, bestModel, DM, exposure, Sangle, Ensemble)
    print(bestModel[3]+np.log(1+np.sum(np.exp(np.array(likelihoods)[:,3]-bestModel[3]))))

def likelihoodPlot_2decays_DMs_better(decay_channels, true_params, radius, num_points, exposure, Sangle,
                                      resolution, load_Ensemble, load_likelihoods, marginalize = True, plotlines = True):
    DM = true_params[0]
    true_decayRate1 = true_params[1]
    true_decayRate2 = true_params[2]
    decayRates1 = np.linspace(0., true_decayRate1 + radius, num_points)
    decayRates2 = np.linspace(0., true_decayRate2 + radius, num_points)
    DM_thresholds = np.zeros(np.size(decay_channels))
    for di in np.arange(np.size(decay_channels)):
        DM_thresholds[di] = np.min(getMasses(decay_channels[di]))
    DMs = getMasses(decay_channels[np.where(DM_thresholds == np.max(DM_thresholds))[0][0]])
    EnsemblePath = Path('Pickles/') / ('Ensemble_'+decay_channels[0]+'_'+decay_channels[1]+'_'+str(true_params[0])+
                                       '_'+str(true_params[1])+'_'+str(true_params[2])+'_'+str(exposure)+'_'+
                                       str(Sangle)+'_'+str(resolution)+'.pk')
    if load_Ensemble:
        try:
            Ensemble = pk.load(open(EnsemblePath, 'rb'))
        except:
            print("Ensemble file doesn't exist, aborting...")
            return
    else:
        Ensemble1 = getEnsemble(decay_channels[0], DM, true_decayRate1, exposure, Sangle)
        Ensemble2 = getEnsemble(decay_channels[1], DM, true_decayRate2, exposure, Sangle)
        Ensemble3 = getBackEnsemble(exposure, Sangle)
        Ensemble = np.concatenate((Ensemble1, Ensemble2, Ensemble3))
        Ensemble = smearEnsemble(Ensemble, resolution)
        print('Ensembles: '+ str(len(Ensemble1))+', '+str(len(Ensemble2))+', '+str(len(Ensemble3)))
        if EnsemblePath.is_file():
            ans = input('Ensemble file already exists. Do you want to overwrite it? (y/n)')
            if ans=='y':
                pk.dump(Ensemble, open(EnsemblePath, 'wb'))
                print('File has been overwritten.')
            else:
                print('File was not overwritten.')
        else:
            pk.dump(Ensemble, open(EnsemblePath, 'wb'))
    likelihoodsPath = Path('Pickles/') / ('loglikelihoods_'+decay_channels[0]+'_'+decay_channels[1]+'_'+
                                          str(true_params[0])+'_'+str(true_params[1])+'_'+str(true_params[2])+'_'+
                                          str(radius)+'_'+str(num_points)+'_'+str(exposure)+'_'+str(Sangle)+'_'+
                                          str(resolution)+'.pk')
    if load_likelihoods:
        try:
            loglikelihoods = pk.load(open(likelihoodsPath, 'rb'))
        except:
            print("Likelihood file doesn't exist, aborting...")
            return
    else:
        loglikelihoods = -np.inf*np.ones((np.size(DMs), np.size(decayRates1), np.size(decayRates2)))
        for m in np.arange(np.size(DMs)):
            E1, PDF1 = getSmearedPDF(decay_channels[0], DMs[m], resolution)
            E2, PDF2 = getSmearedPDF(decay_channels[1], DMs[m], resolution)
            E3, PDF3 = getSmearedBackground(exposure, Sangle, resolution)
            print(f"Mass: {DMs[m]}")
            if True:
                for d1 in np.arange(np.size(decayRates1)):
                    for d2 in np.arange(np.size(decayRates2)):
                        vec = np.array([decayRates1[d1]-true_decayRate1, decayRates2[d2]-true_decayRate2])
                        if np.absolute(np.sum(vec*np.array([1,-1]))/np.sqrt(2)) < 2*radius and np.absolute(np.sum(vec*np.array([1,1]))/np.sqrt(2)) < radius/2:
                            #print('model coordinates: (' + str(DMs[m]) + ',' + str(decayRates1[d1]) + ',' + str(decayRates2[d2]) + ')')
                            Jprob = JPD_2decays(E1, PDF1, E2, PDF2, E3, PDF3, decayRates1[d1], decayRates2[d2], DMs[m], exposure, Sangle, Ensemble)
                            loglikelihoods[m, d1, d2] = Jprob
        if likelihoodsPath.is_file():
            ans = input('Likelihoods file already exists. Do you want to overwrite it? (y/n)')
            if ans=='y':
                pk.dump(loglikelihoods, open(likelihoodsPath, 'wb'))
                print('File has been overwritten.')
            else:
                print('File was not overwritten.')
        else:
            pk.dump(loglikelihoods, open(likelihoodsPath, 'wb'))
    maxJprob = np.max(loglikelihoods)
    indices = np.where(loglikelihoods == maxJprob)
    bestModel = [DMs[indices[0][0]], decayRates1[indices[1][0]], decayRates2[indices[2][0]]]
    trueJprob = JPD_2decays(E1, PDF1, E2, PDF2, E3, PDF3, true_decayRate1, true_decayRate2, DM, exposure, Sangle, Ensemble)
    trueModel = [DM, true_decayRate1, true_decayRate2]
    #bestModel = trueModel
    #maxJprob = trueJprob
    print('Best model:')
    print(bestModel)
    print('True model:')
    print(trueModel)
    sigma5 = 1+scipy.special.erf(-(5)/np.sqrt(2))#5 sigma ellipse
    sigma2 = 1+scipy.special.erf(-(2)/np.sqrt(2))#2 sigma ellipse
    sigmaLikelihoods = sigmaProb(maxJprob, loglikelihoods, deg = 3)
    if marginalize:
        X, Y = np.meshgrid(decayRates1, decayRates2)
        fig, ax = plt.subplots()
        for m in np.arange(np.size(DMs)):
            if True:
                ax.contourf(X, Y, np.transpose(sigmaLikelihoods[m]), [sigma5, 1.], colors = 'purple')
        for m in np.arange(np.size(DMs)):
            if True:
                ax.contourf(X, Y, np.transpose(sigmaLikelihoods[m]), [sigma2, 1.], colors = 'blue')
        '''
        for m in np.arange(np.size(DMs)):
            if DMs[m] == 1000.0:
                sigmaLikelihoods = sigmaProb(maxJprob, loglikelihoods, deg = 2)
                ax.contourf(X, Y, np.transpose(sigmaLikelihoods[m]), [sigma, 1.], colors = 'blue')
        '''
        ax.set_title('Probability of True Model: {:.5%}'.format(sigmaProb(maxJprob, trueJprob))+'\n'
                 + 'True Model: ({:.5} , {:.5}, {:.5})'.format(float(DM), float(true_decayRate1), float(true_decayRate2))+'\n'
                 + 'Best Model: ({:.5} , {:.5}, {:.5})'.format(float(bestModel[0]), float(bestModel[1]), float(bestModel[2])))
        ax.set_ylabel(decay_channels[1]+' Decay Rate [s^-1]')
        ax.set_xlabel(decay_channels[0]+' Decay Rate [s^-1]')

        ax.plot(true_params[1],true_params[2],'kx')

        if plotlines:
            ax.plot(decayRates1, -decayRates1+(true_decayRate1+true_decayRate2+radius/(2*np.sqrt(2))), 'k')
            ax.plot(decayRates1, -decayRates1+(true_decayRate1+true_decayRate2-radius/(2*np.sqrt(2))), 'k')
            ax.plot(decayRates1, decayRates1+(-true_decayRate1+true_decayRate2+2*radius/(np.sqrt(2))), 'k')
            ax.plot(decayRates1, decayRates1+(-true_decayRate1+true_decayRate2-2*radius/(np.sqrt(2))), 'k')
        
        fig.show()

def modelPlot(trueModel, bestModel, DM, exposure, Sangle, Ensemble):
    JfactorDracoDecay = 5.77*(10**24) #(MeV cm^-2 sr^-1)
    secPerYear = 3.154*10**7
#get spectra:
    '''
    E1, PDF1 = getPDF('Pi Eta', DM)
    E2, PDF2 = getPDF('Pi Pi Eta', DM)
    '''
    E1, PDF1 = getSmearedPDF('Pi Eta', DM, .3)
    E2, PDF2 = getSmearedPDF('Pi Pi Eta', DM, .3)
    E3 = np.linspace(startE, 1000, 1000)
    '''
    PDF3 = []
    for energy in E3:
        PDF3.append(background(energy, exposure, Sangle))
    '''
    E3, PDF3 = getSmearedBackground(exposure, Sangle, .3)
#get interpolation functons of decay spectra
    f1 = interp1d(E1, PDF1, kind='linear', fill_value='extrapolate')
    f2 = interp1d(E2, PDF2, kind='linear', fill_value='extrapolate')
#get spectral norms:
    E1 = np.linspace(startE, max(E1), 1000)
    E2 = np.linspace(startE, max(E2), 1000)
    PDF1 = []
    PDF2 = []
    for energy in E1:
        PDF1.append(f1(energy))
    for energy in E2:
        PDF2.append(f2(energy))
    specNorm1 = integrate.simps(PDF1, E1)
    specNorm2 = integrate.simps(PDF2, E2)
    specNorm3 = integrate.simps(PDF3, E3)
#trueModel:
#get expected number of events:
    expectedEvents1 = trueModel[0]/(4*np.pi*DM)*JfactorDracoDecay*specNorm1*exposure*secPerYear*Sangle
    expectedEvents2 = trueModel[1]/(4*np.pi*DM)*JfactorDracoDecay*specNorm2*exposure*secPerYear*Sangle
    expectedEvents3 = specNorm3
    totalExpectedEvents = expectedEvents1 + expectedEvents2 + expectedEvents3
    #totalExpectedEvents = expectedEvents1 + expectedEvents2
#combines background and decays into one PDF based on relative numbers of events
    Emax = max(max(E1), max(E2), max(E3))
    Emin = min(min(E1), min(E2), min(E3))
    if Emin < startE:
        Emin = srartE
    trueE = np.linspace(Emin, Emax, 1000)
    truePDF = []
    for energy in trueE:
        funcValue1 = 0 #necessary because the extrapolated interpolations and background
        funcValue2 = 0 #are unreliable outside domain
        funcValue3 = 0
        if energy >= min(E1) and energy <= max(E1):
            funcValue1 = f1(energy)
        if energy >= min(E2) and energy <= max(E2):
            funcValue2 = f2(energy)
        if energy >= min(E3) and energy <= max(E3):
            funcValue3 = background(energy, exposure, Sangle)
        truePDF.append(((expectedEvents1*funcValue1/specNorm1)
                    +(expectedEvents2*funcValue2/specNorm2)
                    +(expectedEvents3*funcValue3/specNorm3))/totalExpectedEvents)
        #            )/totalExpectedEvents)
    truef = interp1d(trueE, truePDF, kind='linear', fill_value= 0)
    trueNormalization = integrate.simps(truePDF, trueE)
#bestModel:
#get expected number of events:
    expectedEvents1 = bestModel[0]/(4*np.pi*DM)*JfactorDracoDecay*specNorm1*exposure*secPerYear*Sangle
    expectedEvents2 = bestModel[1]/(4*np.pi*DM)*JfactorDracoDecay*specNorm2*exposure*secPerYear*Sangle
    expectedEvents3 = specNorm3
    totalExpectedEvents = expectedEvents1 + expectedEvents2 + expectedEvents3
    #totalExpectedEvents = expectedEvents1 + expectedEvents2
#combines background and decays into one PDF based on relative numbers of events
    Emax = max(max(E1), max(E2), max(E3))
    Emin = min(min(E1), min(E2), min(E3))
    if Emin < startE:
        Emin = srartE
    bestE = np.linspace(Emin, Emax, 1000)
    bestPDF = []
    for energy in bestE:
        funcValue1 = 0 #necessary because the extrapolated interpolations and background
        funcValue2 = 0 #are unreliable outside domain
        funcValue3 = 0
        if energy >= min(E1) and energy <= max(E1):
            funcValue1 = f1(energy)
        if energy >= min(E2) and energy <= max(E2):
            funcValue2 = f2(energy)
        if energy >= min(E3) and energy <= max(E3):
            funcValue3 = background(energy, exposure, Sangle)
        bestPDF.append(((expectedEvents1*funcValue1/specNorm1)
                    +(expectedEvents2*funcValue2/specNorm2)
                    +(expectedEvents3*funcValue3/specNorm3))/totalExpectedEvents)
        #            )/totalExpectedEvents)
    bestf = interp1d(bestE, bestPDF, kind='linear', fill_value= 0)
    bestNormalization = integrate.simps(bestPDF, bestE)
#make plot
    fig, ax = plt.subplots()
    ax.hist(Ensemble, bins = 1000, density=True, label='Photon Histogram')
    ax.plot(bestE, np.array(bestPDF)/bestNormalization, '-', label='Best Model')
    ax.plot(trueE, np.array(truePDF)/trueNormalization, '-', label='True Model')
    ax.set_title('Models and Ensemble')
    ax.set_ylabel('normalized counts/model spectrum')
    ax.set_xlabel('Energy [MeV]')
    ax.legend()
    fig.show()

def sigmaGridPlot_2decays(decayRates1, decayRates2, DMs, exposures, Sangle, resolution, num_ave = 1):
    likelihoods = []
    for decay1 in decayRates1:
        for decay2 in decayRates2:

            if -max(decayRates2)/max(decayRates1)*decay1 + max(decayRates2) <= decay2:
                continue
            if -max(decayRates2)/max(decayRates1)*decay1 + max(decayRates2)/4 >= decay2:
                continue

            print(decay1, decay2)
            Jprob = np.zeros((np.size(DMs), np.size(exposures)))
            Jprob_back = np.zeros((np.size(DMs), np.size(exposures)))
            for Di in np.arange(np.size(DMs)):
                for ei in np.arange(np.size(exposures)):
                    for run in np.arange(num_ave):
                        Ensemble1 = getEnsemble('Pi Eta', DMs[Di], decay1, exposures[ei], Sangle)
                        Ensemble2 = getEnsemble('Pi Pi Eta', DMs[Di], decay2, exposures[ei], Sangle)
                        Ensemble3 = getBackEnsemble(exposures[ei], Sangle)
                        Ensemble = np.concatenate((Ensemble1, Ensemble2, Ensemble3))
                        Ensemble = smearEnsemble(Ensemble, resolution)
                        Jprob[Di, ei] += JPD_2decays(decay1, decay2, DMs[Di], exposures[ei], Sangle, Ensemble)
                        Jprob_back[Di, ei] += JPD_2decays(0., 0., DMs[Di], exposures[ei], Sangle, Ensemble)
                    Jprob[Di, ei] /= num_ave
                    Jprob_back[Di, ei] /= num_ave
            likelihoods.append([decay1, decay2, Jprob, Jprob_back])
    sigma = 1+scipy.special.erf(-(5)/np.sqrt(2))
    epsilon = .010
    fig, ax = plt.subplots()
    x = []
    y = []
    for model in likelihoods:
        isBlue = False
        for Di in np.arange(np.size(DMs)):
            backProb = sigmaProb(model[2][Di, 0], model[3][Di, 0])
            #print(backProb, model[2][Di, 0], model[3][Di, 0])
            '''
            if backProb > sigma and backProb < sigma + epsilon:
                ax.plot(model[0], model[1], 'k.')
            '''
            if backProb > sigma or np.isnan(backProb):
                isBlue = True
                if not model[0] in x:
                    x.append(model[0])
                    y.append(model[1])
                else:
                    y[x.index(model[0])] = np.maximum(y[x.index(model[0])], model[1])
        if isBlue:
            ax.plot(model[0], model[1], 'b.')
        else:
            ax.plot(model[0], model[1], 'k.')
    for decay1 in decayRates1:
        for decay2 in decayRates2:
            if -max(decayRates2)/max(decayRates1)*decay1 + max(decayRates2) <= decay2:
                ax.plot(decay1, decay2, 'k.')
            if -max(decayRates2)/max(decayRates1)*decay1 + max(decayRates2)/4 >= decay2:
                ax.plot(decay1, decay2, 'b.')
    d = np.linspace(np.min(x), np.max(x), 100)
    f = interp1d(x, y, kind='cubic', fill_value= 0)
    ax.plot(d, f(d), 'r-')
    if np.size(DMs) == 1:
        pk.dump([x,f], open(Path('Pickles/') / ('5sigline_'+str(DMs[0])+'DM_'+str(exposures[0])+'exposure.pk'), 'wb'))
    else:
        pk.dump([x,f], open(Path('Pickles/') / ('5sigline_marginalizedDM_'+str(exposures[0])+'exposure.pk'), 'wb'))
    ax.legend(loc='upper right')
    fig.show()

def sigmaLinesPlot_2decays(radii, DMs, exposures, Sangle, resolution, num_ave = 1, partitions = 5):
    discoveryPoints = []
    thetas = np.linspace(0, np.pi/2, partitions + 1)
    for theta in thetas:
        for radius_index in np.flip(np.arange(np.size(radii))):
            print(theta, radii[radius_index])
            decay1 = radii[radius_index]*np.cos(theta)
            decay2 = radii[radius_index]*np.sin(theta)
            Jprob = np.zeros((np.size(DMs), np.size(exposures), np.size(radii)))
            for Di in np.arange(np.size(DMs)):
                for ei in np.arange(np.size(exposures)):
                    for run in np.arange(num_ave):
                        Ensemble1 = getEnsemble('Pi Eta', DMs[Di], decay1, exposures[ei], Sangle)
                        Ensemble2 = getEnsemble('Pi Pi Eta', DMs[Di], decay2, exposures[ei], Sangle)
                        Ensemble3 = getBackEnsemble(exposures[ei], Sangle)
                        Ensemble = np.concatenate((Ensemble1, Ensemble2, Ensemble3))
                        Ensemble = smearEnsemble(Ensemble, resolution)
                        for ri in np.arange(np.size(radii)):
                            Jprob[Di, ei, ri] += JPD_2decays(radii[ri]*np.cos(theta), radii[ri]*np.sin(theta), DMs[Di], exposures[ei], Sangle, Ensemble)
                    Jprob[Di, ei] /= num_ave
            fiveSigmaRadius = 0
            sigma = 1+scipy.special.erf(-(5)/np.sqrt(2))
            for ri in np.flip(np.arange(radius_index)):
                prob = sigmaProb(Jprob[Di, ei, radius_index], Jprob[Di, ei, ri])
                if prob < sigma:
                    fiveSigmaRadius = radii[ri]
                    break
            print(fiveSigmaRadius)
            if radii[radius_index] - 2*(radii[radius_index] - fiveSigmaRadius) < 0:
                discoveryPoints.append([theta, radii[radius_index]])
                break
    fig, ax = plt.subplots()
    for theta in thetas:
        for radius in radii:
            if [theta, radius] in discoveryPoints:
                ax.plot(radius*np.cos(theta), radius*np.sin(theta), 'r.')
            else:
                ax.plot(radius*np.cos(theta), radius*np.sin(theta), 'k.')
    ax.legend(loc='upper right')
    fig.show()

def power(E, P, alpha):
    return P*E**alpha

def backgroundFit(Ensemble, exposure, Sangle):
    fig, ax = plt.subplots()
    hist, bins, patches = ax.hist(Ensemble, bins = 1000, density=False, label='Background Histogram')
    print(hist, np.where(hist != 0., hist, 1.))
    '''
    fit = scipy.stats.linregress(x = np.log(bins[:-1]+(bins[1]-bins[0])/2)[:-700],
                                 y =  np.log(np.where(hist != 0., hist, 0.000001))[:-700])
    '''
    popt, pcov = scipy.optimize.curve_fit(power,
                                          xdata = (bins[:-1]+(bins[1]-bins[0])/2)[:-700],
                                          ydata =  (np.where(hist != 0., hist, 0.000001))[:-700],
                                          sigma = np.sqrt(np.where(hist != 0., hist, 0.000001))[:-700])
    print(popt, pcov)
    perr = np.sqrt(np.diag(pcov))
    ax.plot(bins[:-1], popt[0]*bins[:-1]**popt[1], '-', label='Fit Background')
    ax.plot(bins[:-1], background(bins[:-1], exposure, Sangle)*(bins[1]-bins[0]), '-', label='True Background')
    ax.set_title('Background')
    ax.set_ylabel('Counts')
    ax.set_xlabel('Energy [MeV]')
    ax.legend()
    fig.show()
    print('Fit prefactor : ' + str(popt[0]) + ' +/- ' + str(perr[0]))
    print('True prefactor: ' + str(background(1, exposure, Sangle)*(bins[1]-bins[0])))
    print('Fit exponent: ' + str(popt[1]) + ' +/- ' + str(perr[1]))
    print('True exponent: ' + str(-2))

    

def main():
    '''
    mostLikeleyProb = np.NINF
    mostLikeleyMass = 0
    mostLikeleyDecay = ''
    Ensemble = getEnsemble('Pi Pi Eta', 900.0)
    DMassPi = []
    PiLikelihood = []
    DMassPiPi = 
    PiPiLikelihood = []
    print('Pi Eta')
    for mass in getArchive('Pi Eta'):
        prob = JPDsingle('Pi Eta', mass, 0, Ensemble)
        print(mass, prob)
        if not prob == 0:
            DMassPi.append(mass)
            PiLikelihood.append(prob)
        if not prob == 0 and prob > mostLikeleyProb:
            mostLikeleyProb = prob
            mostLikely = mass
            mostLikeleyDecay = 'Pi Eta'
    print('Pi Pi Eta')
    for mass in getArchive('Pi Pi Eta'):
        prob = JPDsingle('Pi Pi Eta', mass, 0, Ensemble)
        print(mass, prob)
        if not prob == 0:
            DMassPiPi.append(mass)
            PiPiLikelihood.append(prob)
        if not prob == 0 and prob > mostLikeleyProb:
            mostLikeleyProb = prob
            mostLikely = mass
            mostLikeleyDecay = 'Pi Pi Eta'
    print('most likely DM mass: ' + str(mostLikely) + ' via ' + mostLikeleyDecay)

    #print(maximum(DMassPiPi, PiPiLikelihood))
interp1d(DMassPiPi, PiPiLikelihood, kind='linear',
                      fill_value='extrapolate')

    print(maximum(DMassPi, PiLikelihood))
    PiFunc = interp1d(DMassPi, PiLikelihood, kind='linear',
                      fill_value='extrapolate')
    
    archiveMaxPlot(DMassPi, PiLikelihood, PiFunc,
                   DMassPiPi, PiPiLikelihood, PiPiFunc)
    '''
    '''
    DM = 1000.0
    true_decayRate1 = 1.e-26
    true_decayRate2 = 1.e-26
    exposure = 3000
    Sangle = 0.000404322
    resolution = .3
    decayRates1 = np.linspace(0.e-26, 3.e-26, 21)
    decayRates2 = np.linspace(0.e-26, 3.e-26, 21)
    likelihoodPlot_2decays(true_decayRate1, true_decayRate2, decayRates1, decayRates2, DM, exposure, Sangle, resolution)
    '''
    
    decay_channels = ['Pi Pi Eta','K+-']
    true_params = [1050.0, 2.e-26, 2.e-26]
    radius = 5.e-26
    num_points = 300
    exposure = 3000
    Sangle = 0.000404322
    resolution = .3
    load_Ensemble = True
    load_likelihoods = False
    marginalize = True
    plotlines = False
    likelihoodPlot_2decays_DMs_better(decay_channels, true_params, radius, num_points, exposure,
                                      Sangle, resolution, load_Ensemble, load_likelihoods, marginalize, plotlines)
    
    '''
    #DMs = np.linspace(830, 1150, 33)
    DMs = [1000.]
    exposures = [3000.]
    Sangle = 0.000404322
    resolution = .3
    decayRates1 = np.linspace(0.e-26, 0.33e-26, 61)
    decayRates2 = np.linspace(0.e-26, 0.33e-26, 61)
    sigmaGridPlot_2decays(decayRates1, decayRates2, DMs, exposures, Sangle, resolution, num_ave = 30)
    '''
    '''
    DMs = [1000.]
    exposures = [3000.]
    Sangle = 0.000404322
    resolution = .3
    radii = np.linspace(0.e-26, 1.e-26, 21)
    sigmaLinesPlot_2decays(radii, DMs, exposures, Sangle, resolution, num_ave = 1, partitions = 5)
    '''
main()
