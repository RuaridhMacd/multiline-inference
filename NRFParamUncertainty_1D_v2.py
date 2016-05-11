#!/usr/bin/python
# -*- coding: utf-8 -*-

# 1D NRF Parameter Uncertainty Calculator - Python 2.7
# Ruaridh Macdonald, MIT, 2016
# rmacd@mit.edu

# Uses NRFGamma class from:
# StandaloneNRFLineCalculator
# Jayson Vavrek, MIT, 2015
# jvavrek@mit.edu

print "1D NRF Parameter Uncertainty Calculator - Python 2.7"
print "Note: T_Debye still not properly implemented! Doppler-broadened widths will be off."

import sys
import numpy as np
from scipy.stats import norm

import NRFUpdate # Import file with functions to use here
#import NRFPhysics
import NRFmultiLine # Import file with functions to use here
import load_line_data # import file to process resonance and attenuation data
import itertools
#from NRFGamma import NRFGamma # Import NRFGamma class for storing the line properties easily
#import matplotlib.pyplot as plt

# Takes (_z, _a, _Elevel, _Egamma, _Width, _prob, _GSprob, _J0, _Jr, _TDebye, _nDens, _sigmaNRLevel, _sigmaNRGamma, _counter)
# and returns the same as well as: (sigmaInt, sigmaDmax, sigmaNRLevel, sigmaNRGamma, alpha, Delta)

# Timing for debugging and comparisons
import time
startTime = time.time()


# Provide command line option functionality
import argparse
parser = argparse.ArgumentParser(description='Processes a standalone database of NRF gammas and returns lists of branced and neighbouring NRF lines which leak the most information, sorted by attenuation coefficient. \nExample:\
    .......................................................................\
    $./multilineSearch.py -neighE=0.01 -matList=matList.txt -branchOut=5 -neighOut=5 -EMin=1 -EMax=2.5 -bremsMin=2 -bremsMax=9 \
    .......................................................................\
    Produces a lists of the 5 branched and neighbouring pairs with the higher attenuation coefficients, with Elevel between 2 and 9 MeV and EGamma between 1 and 2.5 MeV ')
    
    # runfile('C:/Users/Billy/Dropbox (MIT)/ThesisNRF/multiline-inference-billy/NRFParamUncertainty_1D.py', wdir='C:/Users/Billy/Dropbox (MIT)/ThesisNRF/multiline-inference-billy', args='-neighE=0.01 -matList=matList.txt -branchOut=5 -neighOut=5 -EMin=.2 -EMax=3.5 -bremsMin=.2 -bremsMax=9')
    
# -h and --help options exist by default
parser.add_argument('-neighE', help='Energy gap to qualify as ''neighbours'' (MeV), default = 1KeV ', type=float, default=0.001)
parser.add_argument('-matList', help='File name of list of isotopes to check and their number densities \n Format: A Z numDen*A', type=str)
parser.add_argument('-EMin', help='Detector minimum energy (MeV), default = 0 ', type=float, default=0)
parser.add_argument('-EMax', help='Detector maximum energy (MeV), default = 20 ', type=float, default=20)
parser.add_argument('-bremsMin', help='Photon source minimum energy (MeV), default = 0 ', type=float, default=0)
parser.add_argument('-bremsMax', help='Photon source maximum energy (MeV), default = 20 ', type=float, default=20)
parser.add_argument('-neighOut', help='Lists -neighOut worst neighbouring pairs, ranked by cross section, default = 0 ', type=int)
parser.add_argument('-branchOut', help='Lists -branchOut worst branched pairs, ranked by cross section, default = 0 ', type=int)
parser.add_argument('-plotOn', help='Do you want to plot the results, Yes = 1 ', default = 0, type=int)
args = parser.parse_args()

plotOn = args.plotOn

# Check that material list was given by user
if args.matList == None: sys.exit("User must specify a material file describing the foil and warhead isotopic content")
# Check that definition of 'neighbour' is positive
if args.neighE<= 0 : 
    sys.exit("Neighbour energy gap must be >= 0")
else: neighE = args.neighE
# Check that energy ranges make sense
if args.EMax<=args.EMin or args.EMax<=0: sys.exit("\nMaximum detector energy must be greater than minimum detector energy and positive")
if args.bremsMax<=args.bremsMin or args.bremsMax<=0: sys.exit("\nMaximum Bremsstrahlung energy must be greater than minimum Bremsstrahlung enery and positive")

# ------ ------ ------
# Create list of isotopes to work with and the correspondign number densities
matList, nDensList, thickList = NRFmultiLine.parse_materials(args.matList)
    
# Error checking on material lists (though this will normally be caught during reading)
if len(matList) != len(nDensList): sys.exit("Material input incorrect: different numbers of isotopes and number densities in file")
if len(matList) == 0 : sys.exit("Material input incorrect: must contain at least one material")
    
# Build a matrix of all the non-resonant cross sections up front for use in loops later
NRData = np.loadtxt('nonResonantAttenuation.txt',dtype=float,delimiter='|') # This file has to be carefully formatted
NRData = np.split(NRData,[1],axis=1)
NREnergy = NRData[0] # Vector of energies that the database uses
NRData = NRData[1]   # Non-resonant attenuation data by isotope, z = 1:100
    
# ------ ------ ------
# Import data (read-only) from standalone.dat database of gammas
f = open('standalone.dat','r')
print f
NRfile = open('nonResonantAttenuation.txt','r')
print NRfile

# ------ ------ ------
# Print go statement    
print "\nLooking for significant pairs of lines emitted between %s and %s MeV \nUsing Bremsstrahlung source from %s to %s MeV" %("{:6.3f}".format(args.EMin), "{:6.3f}".format(args.EMax),"{:6.3f}".format(args.bremsMin), "{:6.3f}".format(args.bremsMax))

# ------ ------ ------
# Begin main thread 

# List to contain the NRF lines
isotope_array = load_line_data.load_line_data('mathematicaexport1.csv')
emitList = load_line_data.build_NRFGamma_list(matList, nDensList, thickList, isotope_array, args.EMin, args.EMax, args.bremsMin, args.bremsMax)

# ------ ------ ------
# Calculate counts at detector given the top lines
# Currently we find the highest intensity NRF peak for each isotope and set the smallest one to be 1e4
# This is then used to scale all the other lines

# Find the largest line for each isotope
maxLineCount = [0.0]*len(matList)
maxLineIndx = [0]*len(matList)
        
for i, line in enumerate(emitList):
    if line.counts > maxLineCount[matList.index([line.z,line.a])]:
        maxLineIndx[matList.index([line.z,line.a])] = line.index
        maxLineCount[matList.index([line.z,line.a])] = line.counts
    
  
# Find smallest of these max lines and set source strength so that this line has 1e4 counts
minMaxLineIndex = maxLineCount.index(min(filter(lambda x: x>0,maxLineCount)))

sourceStrength = 1e4/maxLineCount[minMaxLineIndex]
#sourceStrength = 1e10

# Update the count values with the source strength
for line in emitList:
    line.counts = int(line.counts*sourceStrength) #change to always be int. do ceiling #delete emission lines with less than like 100 counts (make this some other variable)
    
# For each line, the number of counts is calculated in the NRFGamma class, assuming Bremms flux = 1 but this is fine for our comparisons (but it makes it trickier to order them)

# ------ ------ ------
# Compare branched and neighbouring pairs vs resonance energy, minimum sigma_NRF and total sigma_NRF
print "Last line index:",emitList[len(emitList)-1].index
branchData, neighData = NRFmultiLine.branchNeigh_compare(emitList,neighE,args.branchOut,args.neighOut,args.plotOn)

# The data lists returned by branchNeigh_compare contain:
#    [0] : The indexes of the branched / neighbouring pairs from emitList, grouped as pairs
#    [1] : The alpha values for the pairs, grouped as pairs
#    [2] : Resonance energies for the lines in the pairs, grouped as pairs for neigh and as a single for branch
#    [3] : Minimum number of counts for the items in a pair
#    [4] : Ordered index of the entries by number of counts, largest -> smallest

# ------ ------ ------
if len(branchData)>0:
# Use numerical methods to estimate properties of interest

    minN   = 0.26
    maxN   = 0.32 #sometimes you need to add a .00000000000000001 or something here to make it work right for the last number.
    deltaN = 0.01
    valsN =  np.arange(minN,maxN+deltaN,deltaN)
    probsNFoil = [0.0]*len(matList) # list of lists of potential number densities for each isotope of interest
     # here each value of Number Density N will be refered to by its index times the delta in this array, while the value is the probability
    for index, isotope in enumerate(matList): 
        probsNFoil[index] = np.divide(np.ones(len(valsN)),len(valsN))
              
    valsN_array = list(list(tup) for tup in itertools.product(valsN,repeat=len(matList))) #returns a list of lists, where each list is one possible combination of N1, N2, N3, ... values
        
    minFoil = 0.8
    maxFoil = 1.2
    deltaFoil = 0.1
    valsFoil = np.arange(minFoil,maxFoil+deltaFoil,deltaFoil)
    probsFoil = np.divide(np.ones(len(valsFoil)),len(valsFoil))
    
    print "Possible values for number density for each isotope:", valsN
    print "Number of pairs:", len(branchData[0])
    
    valsAlpha1Foil = []
    valsAlpha2Foil = []
    probsAlpha1Foil = []
    probsAlpha2Foil = []
    
    branchList = []
    for i in branchData[4]:
        branchList.append(branchData[0][i])
        
    for pairIndex, pair in enumerate(branchList): #iterate over every branched pair. pairIndex is the index of the list of indices of the lines in emitList 
        line1 = emitList[pair[0]] #get reference to NRFGamma class object
        line2 = emitList[pair[1]]    
        
        t1, t2 = NRFUpdate.makeAlphaBranch(line1,line2,valsN_array,matList) # Make all the alpha values, given the combinations of number density guesses
        
        valsAlpha1Foil.append(t1)
        valsAlpha2Foil.append(t2)
        probsAlpha1Foil.append([]) #should use np.product(probsNFoil_array[permIndex]) if probabilities aren't all the same, but that takes a lot longer 
        probsAlpha2Foil.append([])
    
    ####after all of these loops, we now have a list of probabilities and values for both alphas and the foil, named valsFoil, probsAlpha1Foil etc. ###
    
    #iterate over every pair
    #for each branched pair
        #build alpha distributions as function of Ni vectors, set alpha 1 max, min
        #set alpha1 delta
    
    delta = .001
    minCounts = 100
        
    for pairIndex, pair in enumerate(branchList): #iterate over every branched pair. pairIndex is the index of the list of indices of the lines in emitList 
 
        line1 = emitList[pair[0]] #get reference to NRFGamma class object
        line2 = emitList[pair[1]]
        realCounts1 = float(line1.counts)
        realCounts2 = float(line2.counts)
        
        # Break if there aren't enough counts 
        if realCounts1 < minCounts or realCounts2 < minCounts: 
            print "Skipping pair", pairIndex
            continue 
            
        probsAlpha1Foil[pairIndex] = NRFUpdate.makeProbList(probsNFoil)
    
        guessArray = np.zeros([len(valsN)**(len(matList)),len(valsFoil)])
    
        realRatio = realCounts1/realCounts2
        
        # Iterate through each possibility
        for indexAlpha1, guessAlpha1Foil in enumerate(valsAlpha1Foil[pairIndex]):
            indexAlpha2 = indexAlpha1
            guessAlpha2Foil = valsAlpha2Foil[pairIndex][indexAlpha2]   
            
            for indexFoil, guessFoil in enumerate(valsFoil):
                    
                proba1a2x = probsAlpha1Foil[pairIndex][indexAlpha1] * probsFoil[indexFoil]
                      
#                guessRatio = NRFPhysics.guessRatioBranch(guessAlpha1Foil,guessAlpha2Foil,guessFoil,line1.prob,line2.prob)
                guessRatio = (line1.prob*(1.0-np.exp(-guessFoil*guessAlpha1Foil))*guessAlpha2Foil)/(line2.prob*(1.0-np.exp(-guessFoil*guessAlpha2Foil))*guessAlpha1Foil)

                probGuess = (norm.cdf((realRatio+delta-guessRatio)/(guessRatio*(1.0/realCounts1 + 1.0/realCounts2)**.5))-norm.cdf((realRatio-delta-guessRatio)/(guessRatio*(1.0/realCounts1 + 1.0/realCounts2)**.5))) * proba1a2x
            
#                if probGuess > 0.0:
#                    print realRatio, guessRatio, probGuess, realCounts1, realCounts2, emitList[pair[0]].alpha[1], emitList[pair[1]].alpha[1]                
                    
                guessArray[indexAlpha1][indexFoil] = probGuess
        
        # Divide each guess by total sum
#        print "Sum of guess probabilities for this pair:", sumGuesses
        sumGuesses = np.sum(guessArray)        
        if sumGuesses > 0.0:
            guessArray = np.divide(guessArray,sumGuesses)
        
        guessArrayMod = np.reshape(guessArray,np.append((np.repeat(len(probsNFoil[0]),len(matList))),len(probsFoil)))
        
        probsNFoil, probsFoil = NRFUpdate.updateProbBranch(probsNFoil,probsFoil,guessArrayMod,matList)
        

print "\n------ ------ ------ ------ ------ ------\n Branched peaks done \n------ ------ ------ ------ ------ ------\n"
print "Foil number density PDF: "
for i in range(len(matList)):
    print probsNFoil[i]
print "Foil thickness PDF: ", probsFoil
print " \n------ ------ ------ ------ ------ ------ \n"

# ------ ------ ------

print "\n------ ------ ------ ------ ------ ------\n Start neighbouring peaks \n------ ------ ------ ------ ------ ------\n"

if len(neighData) > 0:
            
    minWarhead = 1.8    
    maxWarhead = 2.2
    deltaWarhead = 0.1
    valsWarhead = np.arange(minWarhead,maxWarhead+deltaWarhead,deltaWarhead)
    probsWarhead = np.divide(np.ones(len(valsWarhead)),len(valsWarhead))
        
    minN   = 0.21
    maxN   = 0.39 #sometimes you need to add a .00000000000000001 or something here to make it work right for the last number.
    deltaN = 0.03
    
    delta = 0.001

    valsNWarhead = np.arange(minN,maxN+deltaN,deltaN)
    probsNWarhead = [0.0]*len(matList)
    for index, isotope in enumerate(matList):
        probsNWarhead[index] = np.divide(np.ones(len(valsNWarhead)),len(valsNWarhead))
        
    valsList = []
    for i in range(len(matList)):
        valsList.append(valsN)
    for i in range(len(matList)):
        valsList.append(valsNWarhead)
        
#    valsN_array = list(list(tup) for tup in itertools.product(valsN,repeat=2*len(matList)))
    valsN_array = list(list(tup) for tup in itertools.product(*valsList))
    
    valsAlpha1Foil = []
    valsAlpha2Foil = []
    valsAlpha1Warhead = []
    valsAlpha2Warhead = []
    
    probsAlpha = []
    
    neighList = []
    for i in neighData[4]:
        neighList.append(neighData[0][i])
        
    for pairIndex, pair in enumerate(neighList): #iterate over every branched pair. pairIndex is the index of the list of indices of the lines in emitList 
        line1 = emitList[pair[0]] #get reference to NRFGamma class object
        line2 = emitList[pair[1]]    
        
        t1, t2, t3, t4 = NRFUpdate.makeAlphaNeigh(line1,line2,valsN_array,matList)
        
        probsAlpha.append([])
        valsAlpha1Foil.append(t1)
        valsAlpha2Foil.append(t2)
        
        valsAlpha1Warhead.append(t3)
        valsAlpha2Warhead.append(t4)
    
    for pairIndex, pair in enumerate(neighList): #iterate over every branched pair. pairIndex is the index of the list of indices of the lines in emitList 
#        sumGuesses = 0.0
        line1 = emitList[pair[0]] #get reference to NRFGamma class object
        line2 = emitList[pair[1]]
        realCounts1 = float(line1.counts)
        realCounts2 = float(line2.counts)
        
        # Break if there aren't enough counts 
        if realCounts1 < minCounts or realCounts2 < minCounts: 
            print "Skipping neighbouring pair", pairIndex, "| Counts= ",realCounts1,realCounts2
#            print "Sum of Foil prob = ", np.sum(probsNFoil)/len(matList)
            continue 
        
        probsTemp = []
        for i in range(len(matList)):
            probsTemp.append(probsNFoil[i]) 
        for i in range(len(matList)):
            probsTemp.append(probsNWarhead[i])
#            
        probsAlpha[pairIndex] = NRFUpdate.makeProbList(np.array(probsTemp))
    
        guessArray = np.zeros([len(valsN)**(len(matList)*2),len(valsFoil),len(valsWarhead)])
    
        realRatio = realCounts1/realCounts2
        
        # Iterate through each possibility
        for indexAlpha, guessAlpha1Foil in enumerate(valsAlpha1Foil[pairIndex]):
            guessAlpha2Foil = valsAlpha2Foil[pairIndex][indexAlpha]   
            guessAlpha1Warhead = valsAlpha1Warhead[pairIndex][indexAlpha]
            guessAlpha2Warhead = valsAlpha2Warhead[pairIndex][indexAlpha]
        
            for indexFoil, guessFoil in enumerate(valsFoil):
                
                for indexWarhead, guessWarhead in enumerate(valsWarhead):
                    
                    proba1a2x = probsAlpha[pairIndex][indexAlpha] * probsFoil[indexFoil] * probsWarhead[indexWarhead]
                      
                    guessRatio = line1.sigmaInt[1]/line2.sigmaInt[1] * line1.prob/line2.prob * np.exp(-guessAlpha1Warhead*guessWarhead+guessAlpha2Warhead*guessWarhead) * guessAlpha2Foil/guessAlpha1Foil * (1.0-np.exp(-guessFoil*guessAlpha1Foil))/(1.0-np.exp(-guessFoil*guessAlpha2Foil))
                                           
                    probGuess = (norm.cdf((realRatio+delta-guessRatio)/(guessRatio*(1.0/realCounts1 + 1.0/realCounts2)**.5))-norm.cdf((realRatio-delta-guessRatio)/(guessRatio*(1.0/realCounts1 + 1.0/realCounts2)**.5)))
                        
#                    print guessRatio, realRatio, probGuess, probGuess * proba1a2x

                    guessArray[indexAlpha][indexFoil][indexWarhead] = float(probGuess * proba1a2x)

        # Divide each guess by total sum
#        print "Sum of guess probabilities for this neighbouring pair:", sumGuesses
        sumGuesses = np.sum(guessArray)        
        if sumGuesses > 0.0:
            guessArray = np.divide(guessArray,sumGuesses)
        
        arrayShape = np.append( np.append( np.repeat(len(probsNFoil[0]),2*len(matList)) ,len(probsFoil)) ,len(probsWarhead))
        guessArrayMod = np.reshape(guessArray,arrayShape)
        
        probsNFoil, probsFoil, probsNWarhead, probsWarhead = NRFUpdate.updateProbNeigh(probsNFoil,probsFoil,probsNWarhead,probsWarhead,guessArrayMod,matList)
        
        print "\n ------ ------ ------ ------ ------ ------"
        for i in range(len(matList)):
            print probsNFoil[i]
        print valsN
        
        print "\n",probsFoil
        print valsFoil, "\n"
        
        for i in range(len(matList)):
            print probsNWarhead[i]
        print valsNWarhead
        
        print "\n",probsWarhead
        print valsWarhead
        
        print "Sum of Foil prob = ", np.sum(probsNFoil)/len(matList)

print "\n------ ------ ------ ------ ------ ------\n Neighbouring peaks done \n------ ------ ------ ------ ------ ------\n"

# ------ ------ ------
    
# Close files and tidy up
f.close()
NRfile.close()

print "\n------ ------ ------ ------ ------ ------\n FINAL \n------ ------ ------ ------ ------ ------\n"
for i in range(len(matList)):
    print probsNFoil[i]
print valsN

print "\n",probsFoil
print valsFoil, "\n"

for i in range(len(matList)):
    print probsNWarhead[i]
print valsNWarhead

print "\n",probsWarhead
print valsWarhead
print "\n------ ------ ------ ------ ------ ------\n FINAL \n------ ------ ------ ------ ------ ------\n"
       
# ------ ------ ------
print "Took: %s seconds" %"{:5.3}".format(time.time()-startTime)