# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:04:28 2016

@author: Billy
"""

def load_line_data(fileName):

    import csv 
    isotope_array=[] #will be filled with list of isotopes, which will contain list of data rows, which has dictionary of values
    counter = -1
    prevz = 0
    preva = 0
    
    with open (fileName,'rb') as csvfile: #mathematicaexport1.csv
        reader = csv.DictReader(csvfile,fieldnames=['z','a','e_level','e_gamma','t_half','level_width','brj','br0','j0','jr','sigma_int'])
        for row in reader:
            if row['z'] != prevz or row['a'] != preva:
                prevz = row['z']
                preva = row['a']
                isotope_array.append([])
                counter += 1
                isotope_array[counter].append(row)
            else:
                isotope_array[counter].append(row)
            
    return (isotope_array)        

def build_NRFGamma_list(matList, nDensList, thickList, isotope_array, EMin, EMax, bremsMin, bremsMax): #matList here is just a list of isotopes
    import numpy as np
    from NRFGamma import NRFGamma
    
    gamma_list = []
    TDebye_default = 300
    counter = -1
    
    for isotope_index, isotope in enumerate(isotope_array):
        za = [int(isotope[0]['z']), int(isotope[0]['a'])]
        
        if za in matList:
            #print za
            index = matList.index(za)
            for line in isotope:
                #print float(line['e_gamma'])/1000,line['z'],line['a'],float(line['level_width'])/1000>0 , ~np.isinf(float(line['level_width'])) , float(line['brj'])>0, float(line['e_gamma'])/1000>EMin , float(line['e_gamma'])/1000<EMax , float(line['e_level'])/1000>bremsMin , float(line['e_level'])/1000<bremsMax
                
                if float(line['level_width'])/1000>0 and ~np.isinf(float(line['level_width'])) and float(line['brj'])>0 and float(line['e_gamma'])/1000>EMin and float(line['e_gamma'])/1000<EMax and float(line['e_level'])/1000>bremsMin and float(line['e_level'])/1000<bremsMax: #check conditions for this resonant state
                    counter += 1
                    sigmaNRLevelWarhead = 0
                    sigmaNRGammaWarhead = 0
                    sigmaNRLevelFoil = 0
                    sigmaNRGammaFoil = 0
                    
                    sigmaNRLevelWarheadnoN = []
                    sigmaNRGammaWarheadnoN = []
                    sigmaNRLevelFoilnoN = []
                    sigmaNRGammaFoilnoN = []
                    
                    for i in range(len(matList)) :
                        sigmaNRLevelFoil += get_NR_data(line,matList[i])[0]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at resonance energy
                        sigmaNRGammaFoil += get_NR_data(line,matList[i])[1]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at emitted gamma energy
                        sigmaNRLevelWarhead += get_NR_data(line,matList[i])[0]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at resonance energy
                        sigmaNRGammaWarhead += get_NR_data(line,matList[i])[1]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at emitted gamma energy   
                        #CHANGE TO NOT INCLUDE CONTRIBUTIONS FROM OTHER ISOTOPES EITHER HERE OR IN NRFGamma CLASS
                        
                        sigmaNRLevelFoilnoN.append(get_NR_data(line,matList[i])[0])   #not multiplied by nDensList[i][1] # Non-resonant attenuaton from all isotopes at resonance energy
                        sigmaNRGammaFoilnoN.append(get_NR_data(line,matList[i])[1])    # Non-resonant attenuaton from all isotopes at emitted gamma energy
                        sigmaNRLevelWarheadnoN.append(get_NR_data(line,matList[i])[0]) # Non-resonant attenuaton from all isotopes at resonance energy
                        sigmaNRGammaWarheadnoN.append(get_NR_data(line,matList[i])[1]) # Non-resonant attenuaton from all isotopes at emitted gamma energy 
                        
                    sigmaNRLevel = [sigmaNRLevelWarhead,sigmaNRLevelFoil]
                    sigmaNRGamma = [sigmaNRGammaWarhead,sigmaNRGammaFoil]
                    
                    sigmaNRLevelnoN = [sigmaNRLevelWarheadnoN,sigmaNRLevelFoilnoN]
                    sigmaNRGammanoN = [sigmaNRGammaWarheadnoN,sigmaNRGammaFoilnoN]
                    
                    x = NRFGamma(int(line['z']), int(line['a']), float(line['e_level'])/1000, float(line['e_gamma'])/1000, float(line['level_width'])/1000, float(line['brj']), float(line['br0']), float(line['j0']), float(line['jr']), TDebye_default, nDensList[index], thickList[index], sigmaNRLevel, sigmaNRGamma, counter, _matList=matList, _nDensList=nDensList, _sigmaNRLevelnoN=sigmaNRLevelnoN,_sigmaNRGammanoN=sigmaNRGammanoN)
                        
                    if x.counts > 10**-30:
                        gamma_list.append(x)
                    else:
                        counter -= 1

    #print gamma_list
    return gamma_list
    
def get_NR_data(line, isotope):
    index = line[None].index(str(isotope[0]))
    if line[None][index + 1] != 'NoData' and line[None][index + 2] != 'NoData':
        nrlevel = float(line[None][index + 1])
        nrgamma = float(line[None][index + 2])
        return [nrlevel, nrgamma]
    else: return [0,0]
                    
        
        
        
# old code from main script: 
        
#counter = -1 # Counter for managing indexing in the final list of NRF lines
#
#f.seek(0)
#for line in f: # Running through each line of the NRF database to check if it matches our parameters #this is the section to replace with mathematica data output
#    line = line.split(" ")
#
#    # Grab data and convert strings to ints, floats
#    [z, a] = map(int, [line[0], line[1]]);
#    [Elevel, Egamma, Width, prob, GSprob, J0, Jr, TDebye] = map(float, [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]])
#    
#    # Exclude gammas without valid data, outside the energy ranges or not on the material list
#    if Width>0 and ~np.isinf(Width) and prob>0 and [z,a] in matList and Egamma>args.EMin and Egamma<args.EMax and Elevel>args.bremsMin and Egamma<args.bremsMax:
#        counter += 1 # Add another accepted line to the count
#        
#        nDens = nDensList[matList.index([z,a])]       # Calculate number density * 1e-24 
#        
#        # Search and interpolate for the non-resonant cross section    #not necessary with new data file
#        # Data file has one column of energy labels and then isotopes 1:100
#        levelIndex = NRFmultiLine.find_nearestE(Elevel,NREnergy)# Find energy index of the non-resonant cross sections
#        gammaIndex = NRFmultiLine.find_nearestE(Egamma,NREnergy)    # Find energy index of the non-resonant cross sections
#        
#        sigmaNRLevelWarhead = 0
#        sigmaNRGammaWarhead = 0
#        sigmaNRLevelFoil = 0
#        sigmaNRGammaFoil = 0
#        for i in range(len(matList)) :
#            sigmaNRLevelFoil += NRData[levelIndex,matList[i][0]]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at resonance energy
#            sigmaNRGammaFoil += NRData[gammaIndex,matList[i][0]]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at emitted gamma energy
#            sigmaNRLevelWarhead += NRData[levelIndex,matList[i][0]]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at resonance energy
#            sigmaNRGammaWarhead += NRData[gammaIndex,matList[i][0]]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at emitted gamma energy   
#        sigmaNRLevel = [sigmaNRLevelWarhead,sigmaNRLevelFoil]
#        sigmaNRGamma = [sigmaNRGammaWarhead,sigmaNRGammaFoil]
#        
#        # Build NRFGamma instance
#        x = NRFGamma(z,a,Elevel,Egamma,Width,prob,GSprob,J0,Jr,TDebye,nDens,thickList[matList.index([z,a])],sigmaNRLevel,sigmaNRGamma,counter)          
#        emitList.append(x)


#for index, isotope in enumerate(isotope_array):
#    z = int(isotope[0]['z'])
#    znr = isotope[0][None][index*3]
#    [macl,macg]= [isotope_array[0][0][None][(index*3)+1],isotope_array[0][0][None][(index*3)+2]]
#    print z,znr,macl,macg