#!/usr/bin/python
# -*- coding: utf-8 -*-

# NRFGamma class from StandaloneNRFLineCalculator
# Jayson Vavrek, MIT, 2015
# jvavrek@mit.edu

# Modified by Ruaridh Macdonald
# rmacd@mit.edu

import math

# Build quick class structure
class NRFGamma:
    def __init__(self, _z, _a, _Elevel, _Egamma, _Width, _prob, _GSprob, _J0, _Jr, _TDebye, _nDens, _thickness, _sigmaNRLevel, _sigmaNRGamma, _counter, _matList=[], _nDensList=[], _sigmaNRGammanoN=[], _sigmaNRLevelnoN=[]):
        self.z      = _z
        self.a      = _a
        self.Elevel = _Elevel  # Energy of the resonant level
        self.Egamma = _Egamma  # Energy of the emitted gamma
        self.Width  = _Width   # Width Gamma_r of the resonant level
        self.prob   = _prob    # Branching ratio brj of decay from resonance to final state
        self.GSprob = _GSprob  # Branching ratio b0r of decay from resonance to ground sate
        self.J0     = _J0      
        self.Jr     = _Jr      
        self.TDebye = _TDebye
        self.nDens  = _nDens   # atom / cm^2 * 1e-24
        self.thickness = _thickness # [0]=warhead thickness, [1]=foil thickness
        self.index = _counter
        self.matList = _matList
        self.nDensList = _nDensList
        self.sigmaNRLevel = _sigmaNRLevel
        self.sigmaNRGamma = _sigmaNRGamma
        self.sigmaNRGammanoN = _sigmaNRGammanoN #list of warhead and foil lists of NR gamma contributions not yet multiplied by number densities of each isotope in matList
        self.sigmaNRLevelnoN = _sigmaNRLevelnoN

        # Calculate the energy-integrated cross section
        g = (2.0*self.Jr+1)/(2.0*(2.0*self.J0+1))
        hbarc = 197.327e-15 # MeV m
        self.sigmaInt = 1.0e34 * 2.0 * (math.pi)**2 * g * (hbarc/self.Elevel)**2 * self.Width * self.prob * self.GSprob # eV b
        #print self.sigmaInt, self.Elevel, self.z,self.a
        self.sigmaInt = [self.sigmaInt/self.prob * self.nDens[0] , self.sigmaInt/self.prob * self.nDens[1]]
        
        # and the Doppler-broadened peak height
        Mc2 = self.a * 931.454 # MeV
        kB = 8.6173e-11 # MeV/K
        self.Delta = self.Elevel * math.sqrt(2*kB*300.0/Mc2) # FIXME this should use Teff
        self.sigmaDmax = 1.0e28 * 2.0 * (math.pi)**(3.0/2.0) * g * (hbarc/self.Elevel)**2 * self.prob * self.GSprob * self.Width / self.Delta # b
        
        # Calculate the alpha factor : mu_NRF(Elevel) + mu_NR(Elevel) + 2*mu_NR(Egamma) [Warhead,Foil]
        #self.sigmaNRLevel = [ _sigmaNRLevel[0] *self.Delta *self.nDens[0] , _sigmaNRLevel[1] *self.Delta *self.nDens[1] ] # Multiply non-resonant cross sections by resonance width [Warhead,Foil]
        #self.sigmaNRGamma = [ _sigmaNRGamma[0] *self.Delta *self.nDens[0] , _sigmaNRGamma[1] *self.Delta *self.nDens[1] ] # Note that they already have number density incorporated [Warhead,Foil]
       
       #self.alpha = [self.sigmaInt[0] + math.log(1/self.prob) + self.sigmaNRLevel[0] , self.sigmaInt[1] + math.log(1/self.prob) + self.sigmaNRLevel[1] + 2*self.sigmaNRGamma[1]]
        self.alpha = [self.sigmaInt[0] + self.sigmaNRLevel[0] , self.sigmaInt[1] + self.sigmaNRLevel[1] + 2*self.sigmaNRGamma[1]]
            
        self.counts = math.exp(-self.alpha[0]*self.thickness[0]) * self.prob * self.sigmaInt[1] / self.alpha[1] * (1 - math.exp(-self.alpha[1]*self.thickness[1]))