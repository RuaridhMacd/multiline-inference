# -*- coding: utf-8 -*-
"""
Created on Wed May 04 15:39:12 2016

@author: Ruaridh
"""

import numpy as np
import itertools

# ------ ------ ------

def alphaBranch(line,valsN_array,matList):
    
    valsAlpha = []
    
    sigmaNRF = line.sigmaInt[1]/line.nDens[1]
    sigmaNR_levels = line.sigmaNRLevelnoN[1]
    sigmaNR_gammas = line.sigmaNRGammanoN[1]
    
    for permIndex, permutationOfNs in enumerate(valsN_array): #each permutationN is a list of N values for each isotope in matList
        sigmaNR_level = 0.0
        sigmaNR_gamma = 0.0
        nDensIsotope  = permutationOfNs[matList.index([line.z,line.a])]
        
        for matIndex, mat in enumerate(matList):
            sigmaNR_level += sigmaNR_levels[matIndex]*permutationOfNs[matIndex] #multiply cross section times number density 
            sigmaNR_gamma += sigmaNR_gammas[matIndex]*permutationOfNs[matIndex]
                
        valsAlpha.append(sigmaNRF*nDensIsotope+sigmaNR_level+2*sigmaNR_gamma) #estimated alpha values for this selection of nDens possibilities, 

    return valsAlpha
    
# ------ ------ ------

def alphaNeigh(line,valsN_array,matList):
    
    valsAlpha = []
    
    sigmaNRF = line.sigmaInt[0]/line.nDens[0]
    sigmaNR_levels = line.sigmaNRLevelnoN[0]
    
    for permIndex, permutationOfNs in enumerate(valsN_array): #each permutationN is a list of N values for each isotope in matList
        sigmaNR_level = 0.0
        nDensIsotope  = permutationOfNs[matList.index([line.z,line.a])]
        
        for matIndex, mat in enumerate(matList):
            sigmaNR_level += sigmaNR_levels[matIndex]*permutationOfNs[matIndex+len(matList)] #multiply cross section times number density 
                
        valsAlpha.append(sigmaNRF*nDensIsotope+sigmaNR_level) #estimated alpha values for this selection of nDens possibilities, 

    return valsAlpha
    
# ------ ------ ------

def makeAlphaBranch(line1,line2,valsN_array,matList):
    
    valsAlpha1 = alphaBranch(line1,valsN_array,matList)
    valsAlpha2 = alphaBranch(line2,valsN_array,matList)
    
    return valsAlpha1, valsAlpha2

# ------ ------ ------

def makeAlphaNeigh(line1,line2,valsN_array,matList):
    
    valsAlpha1Foil = alphaBranch(line1,valsN_array,matList)
    valsAlpha2Foil = alphaBranch(line2,valsN_array,matList)
    
    valsAlpha1Warhead = alphaNeigh(line1,valsN_array,matList)
    valsAlpha2Warhead = alphaNeigh(line2,valsN_array,matList)
        
    return valsAlpha1Foil, valsAlpha2Foil, valsAlpha1Warhead, valsAlpha2Warhead
    
# ------ ------ ------

def makeProbList(probsN):
    
    tupList = np.array(list(tup for tup in itertools.product(*probsN))) # Produce an np array of all the combinations of variables
    tupProb = np.prod(tupList,axis=1)                                   # Take product of each tuple = prob[combination]
    
    return tupProb
 
# ------ ------ ------

def updateProbBranch(probsNFoil,probsFoil,guessArrayMod,matList):
    
    lowerThreshold = 1e-20 
    
    probsNFoil.append(probsFoil) # Append the foil thickness PDFs to make it all one matrix

    counter = -1
    for i in itertools.combinations(reversed(range(len(guessArrayMod.shape))),len(matList)):
        counter += 1
        probsNFoil[counter] = np.sum(guessArrayMod,axis=i) # For each PDF, sum across the other axes to get the new PDF
    probsFoil = probsNFoil.pop()
    
    for j in range(len(matList)): 
        probsNFoil[j][probsNFoil[j]<lowerThreshold] = 0.0 # Zero the values which are very unlikely 
    probsFoil[probsFoil<lowerThreshold] = 0.0
    
    return probsNFoil, probsFoil

# ------ ------ ------
        
def updateProbNeigh(probsNFoil,probsFoil,probsNWarhead,probsWarhead,guessArrayMod,matList):
    
    lowerThreshold = 1e-20
    
    if np.sum(guessArrayMod)>0.0:                               # Only update if we won't divide by zero
        probsTemp = []
        for i in range(len(matList)):
            probsTemp.append(probsNFoil[i])                     # Combine the foil and warhead numDens PDFs
        for i in range(len(matList)):
            probsTemp.append(probsNWarhead[i])
        probsTemp.append(probsFoil)                             # Append the thickness PDFs to make it all one matrix
        probsTemp.append(probsWarhead)
        
        counter = -1
        for i in itertools.combinations(reversed(range(len(guessArrayMod.shape))),len(matList)*2+1):
            counter += 1
            probsTemp[counter] = np.sum(guessArrayMod,axis=i)   # For each PDF, sum across the other axes to get the new PDF
            
        probsWarhead = probsTemp.pop()                          # Separate out the thickness PDFs in reverse order
        probsFoil = probsTemp.pop()
        probsNWarhead = probsTemp[len(matList):2*len(matList)]  # Assign the number density PDFs
        probsNFoil = probsTemp[0:len(matList)]
        
        for j in range(len(matList)):
            probsNFoil[j][probsNFoil[j] < lowerThreshold] = 0.0        # Zero the values which are very unlikely 
            probsNWarhead[j][probsNWarhead[j] < lowerThreshold] = 0.0
        probsFoil[probsFoil < lowerThreshold] = 0.0
        probsWarhead[probsWarhead < lowerThreshold] = 0.0
        
    else:
        print "\nToo small to update"

    return probsNFoil, probsFoil, probsNWarhead, probsWarhead