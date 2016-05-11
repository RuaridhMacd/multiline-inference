#!/usr/bin/python
# -*- coding: utf-8 -*-

# Functions useful for NRF multiline analysis - Python 2.7
# Ruaridh Macdonald, MIT, 2016
# rmacd@mit.edu

import numpy as np
import matplotlib.pyplot as plt
import itertools
import math

import scipy.stats as stat

# ------ ------ ------

# Function to parse material input and give list of isotopes and their densities in the foil and warhead
def parse_materials(fileName):
    if type(fileName) != str:
        print('parse_material ERROR: \nMaterial list name must be a string')
        return
        
    matList = []   # List of atomic and mass numbers for materials in object
    nDensList = [] # (atom / cm**3) * A * 1e-24 
    thickList = [] # [warhead,foil] thicknesses (cm)
    
    matFile = open(fileName,'r')
    print matFile
    
    for isotope in matFile:
        isotope = isotope.strip()
        isoColumn = isotope.split(' ')
        
        if isoColumn[0]=="#": continue # Ignore comment lines in the material description
        
        matList.append([int(isoColumn[0]),int(isoColumn[1])])
        nDensList.append([float(isoColumn[2]),float(isoColumn[3])])
        thickList.append([float(isoColumn[4]),float(isoColumn[5])])
        
    matFile.close()
    
    return (matList,nDensList,thickList)
    
# ------ ------ ------
    
# Function to do energy searches
def find_nearestE(ELevel,NREnergy):
    idx = (np.abs(NREnergy-ELevel)).argmin()
    return idx

# ------ ------ ------

# Function to iterate through pairs in emitList and find the most significant neighbouring and branched peaks, given their NRF cross sections
def branchNeigh_compare(emitList,deltaNeigh,numBranch,numNeigh,plotOn):
    
    indexListNeigh = [] 
    alphaListNeigh = []
    ratioListNeigh = []
    levelListNeigh = []
    medsigIntListNeigh = []
    minCListNeigh = []
    cListPlotNeigh =  []

    indexListBranch = []
    eListPlotBranch = []
    sListPlotBranch = []
    sTListPlotBranch = []
    aListPlotBranch = []
    rListTempBranch = []
    cListPlotBranch = []
    branchData = []
    
    # Use itertools.combinations to make the smallest list of all combinations of pairs
    # It doesn't include any duplicates and doesn't pair items with themselves
    for pair in itertools.combinations(emitList,2):
        # Neighbouring pairs
        if abs(pair[0].Elevel-pair[1].Elevel)<=deltaNeigh and pair[0].Elevel != pair[1].Elevel :
            
            indexListNeigh.append([pair[0].index, pair[1].index]) # Store indexes of neighbouring pair
            
            alphaListNeigh.append(min([pair[0].alpha[1] , pair[1].alpha[1]]))   # Alpha values of pair
            ratioListNeigh.append(pair[0].alpha[1]/pair[1].alpha[1])            # Alpha ratios of pair
            levelListNeigh.append(np.median([pair[0].Elevel , pair[1].Elevel])) # Resonance energies of pair
            medsigIntListNeigh.append(np.median([pair[0].sigmaInt[1]+math.log(1/pair[0].prob), pair[1].sigmaInt[1]+math.log(1/pair[1].prob)]))  # Median NRF cross section of pair
            minCListNeigh.append(np.min([pair[0].sigmaInt[1]+math.log(1/pair[0].prob) , pair[1].sigmaInt[1]+math.log(1/pair[1].prob)]))    # Minimum NRF cross section of pair
            cListPlotNeigh.append(np.min([pair[0].counts,pair[1].counts]))
        # Branched pairs
        elif pair[0].Elevel == pair[1].Elevel and pair[0].a==pair[1].a and pair[0].z==pair[1].z  :
            indexListBranch.append([pair[0].index , pair[1].index])
            aListPlotBranch.append(max([pair[0].alpha[1] , pair[1].alpha[1]]))
            rListTempBranch.append(pair[0].alpha[1]/pair[1].alpha[1])
            eListPlotBranch.append(pair[0].Elevel)
            sTListPlotBranch.append(pair[0].sigmaInt[1]+math.log(1/pair[0].prob))
            sListPlotBranch.append(min([pair[0].sigmaInt[1] , pair[1].sigmaInt[1]]))
            cListPlotBranch.append(np.min([pair[0].counts,pair[1].counts]))
    
    if len(rListTempBranch) == 0:#
        print('\nNo branched pairs')
        if plotOn == 1:
            plt.figure(1)
            plt.clf()
            plt.suptitle('No branched pairs', fontsize=24)
        branchData = []
        
    else:
        # Find and print list of most significant branched lines based on counts
        if numBranch != None:
            
            minCListBranch = np.array(cListPlotBranch)
            minCListBranchIndex = np.argsort(minCListBranch)[::-1]
            
            branchData = [indexListBranch,aListPlotBranch,eListPlotBranch,minCListBranch,minCListBranchIndex]

            if len(minCListBranch) < numBranch: 
                print('\nOnly %i branched pairs found' %len(minCListBranch))
                print('The %i branched pairs with largest total NRF attenuation coefficient: \n  Isotope    mu_NRF[b] branch ratio Elevel[MeV]    EGamma[MeV]          alpha     alpha_ratio  counts' % numBranch)
                for i in range(len(minCListBranch)):
                    print "[%s , %s] %s  [%s , %s] %s    [%s , %s] [%s , %s] %s %s" % ("{:3.0f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].z), "{:3.0f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].a), "{:6.3f}".format(minCListBranch[minCListBranchIndex[i]]), "{:4.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].prob), "{:4.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].prob), "{:8.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].Elevel), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].Egamma), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].Egamma), "{:6.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].alpha[1]), "{:6.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].alpha[1]), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].alpha[1]/emitList[indexListBranch[minCListBranchIndex[i]][1]].alpha[1]), "{:10.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].counts))
            else:
                print('\nFound %i branched pairs' %len(minCListBranch) )
                print('The %i branched pairs with largest total NRF attenuation coefficient: \n  Isotope    mu_NRF[b] branch ratio Elevel[MeV]    EGamma[MeV]          alpha     alpha_ratio  counts' % numBranch)
                for i in range(numBranch):
                    print "[%s , %s] %s  [%s , %s] %s    [%s , %s] [%s , %s] %s %s" % ("{:3.0f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].z), "{:3.0f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].a), "{:6.3f}".format(minCListBranch[minCListBranchIndex[i]]), "{:4.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].prob), "{:4.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].prob), "{:8.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].Elevel), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].Egamma), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].Egamma), "{:6.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].alpha[1]), "{:6.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].alpha[1]), "{:6.3f}".format(emitList[indexListBranch[minCListBranchIndex[i]][0]].alpha[1]/emitList[indexListBranch[minCListBranchIndex[i]][1]].alpha[1]), "{:10.2f}".format(emitList[indexListBranch[minCListBranchIndex[i]][1]].counts))
              
                if plotOn == 1:
                    plt.figure(1)
                    plt.clf()
                    plt.suptitle('Branched emissions', fontsize=18)
                    plt.subplot(221)
                    plt.semilogy(eListPlotBranch,aListPlotBranch,'.k')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Resonance energy (MeV)',fontsize=18)
                    plt.ylabel('$\\alpha_{0,1}$ (b)',fontsize=24)
                    plt.subplot(222)
                    plt.semilogy(eListPlotBranch,rListTempBranch,'.k')
                    plt.semilogy([0,max(eListPlotBranch)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Resonance energy (MeV)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{0,2}}$',fontsize=40)
                    plt.subplot(223)
                    plt.semilogy(sListPlotBranch,rListTempBranch,'.k')
                    plt.semilogy([0,max(sListPlotBranch)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Min B$*\sigma_{NRF}$ (b)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{0,2}}$',fontsize=40)
                    plt.subplot(224)
                    plt.semilogy(minCListBranch,rListTempBranch,'.k')
                    plt.semilogy([0,max(minCListBranch)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('$\sigma_{NRF}$ (b)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{0,2}}$',fontsize=40)
                    plt.show()
            
    if len(ratioListNeigh) == 0:
        print('\nNo neighbouring pairs')
        neighData = []
        if plotOn == 1:
            plt.figure(2)
            plt.clf()
            plt.suptitle('No neighbouring pairs', fontsize=24)
    else:
        # Print list of most significant neighbouring lines, by counts
        if numNeigh != None:
            
            minCListNeigh = np.array(cListPlotNeigh)
            minCListNeighIndex = np.argsort(minCListNeigh)[::-1]

            neighData = [indexListNeigh,alphaListNeigh,levelListNeigh,minCListNeigh,minCListNeighIndex]

            if len(minCListNeigh) < numNeigh: 
                print('\nOnly %i neighbouring pairs found' %len(minCListNeigh))
                print('The %i neighbouring pairs with largest minimum( mu_NRF of pair ): \n  Istope_1   Isotope_2  mu_NRF[b]  ELevel [MeV]    EGamma [MeV]        alpha     alpha_ratio counts' % numNeigh)
                for i in range(len(minCListNeigh)):
                    print "[%s , %s] [%s , %s] %s  [%s , %s] [%s , %s] [%s , %s] %s %s" %("{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].z),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].a),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].z),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].a), "{:6.4f}".format(minCListNeigh[minCListNeighIndex[i]]), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].Elevel), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].Elevel), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].Egamma), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].Egamma), "{:6.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].alpha[1]), "{:6.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].alpha[1]), "{:6.3f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].alpha[1]/emitList[indexListNeigh[minCListNeighIndex[i]][1]].alpha[1]), "{:10.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].counts)  )
            else:
                print('\nFound %i neighbouring pairs' % len(minCListNeigh))
                print('The %i neighbouring pairs with largest minimum( mu_NRF of pair ): \n  Istope_1   Isotope_2  mu_NRF[b]  ELevel [MeV]    EGamma [MeV]        alpha     alpha_ratio counts' % numNeigh)
                for i in range(numNeigh):
                    print "[%s , %s] [%s , %s] %s  [%s , %s] [%s , %s] [%s , %s] %s %s" %("{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].z),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].a),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].z),"{:3.0f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].a), "{:6.4f}".format(minCListNeigh[minCListNeighIndex[i]]), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].Elevel), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].Elevel), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].Egamma), "{:5.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].Egamma), "{:6.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].alpha[1]), "{:6.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].alpha[1]), "{:6.3f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][0]].alpha[1]/emitList[indexListNeigh[minCListNeighIndex[i]][1]].alpha[1]), "{:10.2f}".format(emitList[indexListNeigh[minCListNeighIndex[i]][1]].counts) )
               
                if plotOn == 1:
                    plt.figure(2)
                    plt.clf()
                    plt.suptitle('Neighbouring emissions', fontsize=18)
                    plt.subplot(221)
                    plt.semilogy(levelListNeigh,alphaListNeigh,'.k')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Resonance Energy (MeV)',fontsize=18)
                    plt.ylabel('Min$(\\alpha_{0,1},\\alpha_{2,3})$',fontsize=18)
                    plt.subplot(222)
                    plt.semilogy(levelListNeigh,ratioListNeigh,'.k')
                    plt.semilogy([0,max(levelListNeigh)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Resonance Energy (MeV)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{2,3}}$',fontsize=40)
                    plt.subplot(223)
                    plt.semilogy(medsigIntListNeigh,ratioListNeigh,'.k')
                    plt.semilogy([0,max(medsigIntListNeigh)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Median $\sigma_{NRF}$ (b)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{2,3}}$',fontsize=40)
                    plt.subplot(224)
                    plt.semilogy(minCListNeigh,ratioListNeigh,'.k')
                    plt.semilogy([0,max(minCListNeigh)],[1,1],'--r')
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    plt.xlabel('Min $\sigma_{NRF}$ (b)',fontsize=18)
                    plt.ylabel('$\\frac{\\alpha_{0,1}}{\\alpha_{2,3}}$',fontsize=40)
                    plt.show()

    return branchData,neighData