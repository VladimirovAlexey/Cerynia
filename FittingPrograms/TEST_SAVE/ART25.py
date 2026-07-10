#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""
#######################################
# importing libraries
#######################################
import os
THIS_DIR= os.path.dirname(__file__)
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"
HARPY_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/harpy/"

MODEL_DIR_ART25 = ATMDE_DIR+"Models/ART25/Replica-files/"

import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

#%%
import Cerynia
import harpy

path_to_constants=THIS_DIR+"/INI/ART25_main.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5004, 0.073018, 0.038048, 0.0]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]
inARRAY_FF=[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)
harpy.setNPparameters_uTMDFF(inARRAY_FF)

#%%
rSet=Cerynia.aTMDeReplicaSet.from_json(MODEL_DIR_ART25+"ART25_main.json")

rSet.set(0)
#%%
### read the list of files and return the list of DataSets
def loadThisDataSIDIS(listOfNames):    
    
    path_to_data=ROOT_DIR+"DataLib/SIDIS/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=Cerynia.DataSet.from_csv(path_to_data+name+".csv")        
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    
    path_to_data=ROOT_DIR+"DataLib/DY/"
    path_to_dataW=ROOT_DIR+"DataLib/DY_W/"
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=Cerynia.DataSet.from_csv(path_to_dataW+name+".csv")        
        elif("_A4" in name):
            loadedData=Cerynia.DataSet.from_csv(path_to_dataA+name+".csv")   
        else:
            loadedData=Cerynia.DataSet.from_csv(path_to_data+name+".csv")        
        dataCollection.append(loadedData)   
        

    return dataCollection

#%%
##################Cut function
def cutFuncSIDIS(p):          
    if p["z_avg"]>0.8:
        return False , p
    ## bins with low z drop
    if p["z_avg"]<0.2:
        return False , p
    
    par=1.0
    if p["xSec"]<0.00000001:
        err=1
        delta=1
    else:
        ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
        err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec    
        gamma2=(2.0*p["M_target"]*p["x_avg"]/p["Q_avg"])**2
        rho2=(p["M_product"]/p["z_avg"]/(p["Q_avg"]))**2
        qT=p["pT_avg"]/p["z_avg"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
        delta=qT/(p["Q_avg"])
        
        ### compute the largest possible qT (approximate)
        gamma2WORST=(2.0*p["M_target"]*p["x_max"]/p["Q_avg"])**2
        # it is definitely not a TMD point
        if gamma2WORST*rho2>1:
            return False , p
        qTWORST=p["pT_max"]/p["z_min"]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))

        ## drop if qT>Q/2
        if qTWORST>p["Q_avg"]/2:
            return False , p

    ### drop Q<2
    if p["Q_avg"]<2 :
        return False , p
    
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%    
def cutFuncDY(p):
    par=0.5
    
    if(p["xSec"]>0):
        err=Cerynia.Point.uncorrTotal(p)/p["xSec"]
    else:
        err=100.
    delta=p["qT_avg"]/p["Q_avg"]
    
    if(p["id"][0] == "E"):
        delta=p["qT_avg"]/p["Q_max"]
    
    if("run1-W" in p["id"]):
        delta=p["qT_min"]/(p["Q_min"]+5.)
    
    
    if(p["id"][0:4] == "E605"):
        if(p["Q_min"]==10.5):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E772"):
        if(p["Q_min"]<10):#these bins seems broken
            return False , p
    elif(p["id"][0:4] == "E615"):
        if(9<p["Q_avg"]<11.2):#UPSILON resonance-bin
            return False , p
    elif(p["id"][0:4] == "E228"):
        if(9<p["Q_avg"]<11):#UPSILON resonance-bin
            return False , p
    else:
        if(9<p["Q_avg"]<11):#UPSILON resonance-bin
            return False , p    
    
    return ((delta<0.25 and p["qT_avg"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    
   

#%%
### Loading the SIDIS data set
theData=Cerynia.DataMultiSet(loadThisDataSIDIS([
                      'hermes3D.p.pi+','hermes3D.p.pi-',
                      'hermes3D.d.pi+','hermes3D.d.pi-',
                      'hermes3D.p.k+','hermes3D.p.k-',
                      'hermes3D.d.k+','hermes3D.d.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.cut(cutFuncSIDIS) 
setSIDIS.prepare()

print(setSIDIS)

#%%
### Loading the DY data set
theData=Cerynia.DataMultiSet(loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy', 
                          'PHE200', 'STAR510', 
                          'E288-200', 'E288-300', 'E288-400', 
                          'E772',
                          'E605',
                          'D01_W','CDF1_W'
                          ]))

setDY=theData.cut(cutFuncDY) 
setDY.prepare()

print(setDY)

#%%
Cerynia.harpyInterface.print_chi2_table(setDY,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setSIDIS,sys_shift=False)