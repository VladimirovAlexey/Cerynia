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
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet
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
import DataProcessor.ArtemideReplicaSet

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR_ART25+"ART25_main.rep")
    
rSet.SetReplica(0)
#%%
### read the list of files and return the list of DataSets
def loadThisDataSIDIS(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolSIDIS/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection

#%%
### read the list of files and return the list of DataSets
def loadThisDataDY(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=ROOT_DIR+"DataLib/unpolDY/"
    path_to_dataW=ROOT_DIR+"DataLib/unpolW/"
    path_to_dataA=ROOT_DIR+"DataLib/DY_angular/"
    
    
    dataCollection=[]
    for name in listOfNames:
        if(name[-1]=="W"):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataW+name+".csv")
        elif("_A4" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        elif("_Auu" in name):
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_dataA+name+".csv")
        else:
            loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        

    return dataCollection

#%%
##################Cut function
def cutFunc(p):
    
    if p["type"]=="SIDIS":        
        if p["<z>"]>0.8:
            return False , p
        ## bins with low z drop
        if p["<z>"]<0.2:
            return False , p
        
        par=1.0
        if p["xSec"]<0.00000001:
            err=1
            delta=1
        else:
            ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
            err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec    
            gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
            rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
            qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
            delta=qT/(p["<Q>"])
            
            ### compute the largest possible qT (approximate)
            gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
            # it is definitely not a TMD point
            if gamma2WORST*rho2>1:
                return False , p
            qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))
    
            ## drop if qT>Q/2
            if qTWORST>p["<Q>"]/2:
                return False , p
    
        ### drop Q<2
        if p["<Q>"]<2 :
            return False , p
        
        #### comensate new definitions
        p["thFactor"]=p["thFactor"]*(p["z"][1]-p["z"][0])*\
                (p["x"][1]-p["x"][0])*(p["pT"][1]-p["pT"][0])*(p["Q"][1]-p["Q"][0])
        
        #### checking new name convention
        if(p["process"]==[1,1,1,2002]): p["process"]=[1,12,1,2001]
        if(p["process"]==[1,1,-1,2002]): p["process"]=[1,12,-1,2001]
        if(p["process"]==[1,1,2,2002]): p["process"]=[1,12,2,2001]
        if(p["process"]==[1,1,-2,2002]): p["process"]=[1,12,-2,2001]
        if(p["process"]==[1,1,1,2103]): p["process"]=[1,12,12,2001]
        if(p["process"]==[1,1,-1,2103]): p["process"]=[1,12,-12,2001]
        
        return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p
    elif p["type"]=="DY":
        par=0.5
        
        #  for artemide v3.    
        # p["process"]=[p["process"][0],p["process"][2],1,1]
        if(len(p["process"])==3):        
                print("UNKNOWN PROCESS IN ARTEMIDE 3"+str(p["process"]))
        
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
        delta=p["<qT>"]/p["<Q>"]
        
        if(p["id"][0] == "E"):
            delta=p["<qT>"]/p["Q"][1] 
        
        if("run1-W" in p["id"]):
            delta=p["qT"][0]/(p["Q"][0]+5.)
        
        
        if(p["id"][0:4] == "E605"):
            if(p["Q"][0]==10.5):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E772"):
            if(p["Q"][0]<10):#these bins seems broken
                return False , p
        elif(p["id"][0:4] == "E615"):
            if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E228"):
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
        else:
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p    
        
        #### comensate new definitions
        p["thFactor"]=p["thFactor"]*(p["y"][1]-p["y"][0])*(p["qT"][1]-p["qT"][0])*(p["Q"][1]-p["Q"][0])
        
        #### checking new name convention
        if(p["process"]==[2,1,1,101]): p["process"]=[2,1,64029,1]

        
        return ((delta<0.25 and p["<qT>"]<10.) or (delta<0.25 and par/err*delta**2<1)) , p
    
   

#%%
### Loading the SIDIS data set
theData=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSIDIS([
                      'hermes.p.vmsub.zxpt.pi+','hermes.p.vmsub.zxpt.pi-',
                      'hermes.d.vmsub.zxpt.pi+','hermes.d.vmsub.zxpt.pi-',
                      'hermes.p.vmsub.zxpt.k+','hermes.p.vmsub.zxpt.k-',
                      'hermes.d.vmsub.zxpt.k+','hermes.d.vmsub.zxpt.k-',
                      'compass.d.h+','compass.d.h-']))

setSIDIS=theData.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setSIDIS.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSIDIS.sets]), 'points.')

#%%
### Loading the DY data set
theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataDY([
                          'CDF1', 'CDF2', 'D01', 'D02', 'D02m', 
                          #'A7-00y10', 'A7-10y20','A7-20y24', 
                          'A8-00y04', 'A8-04y08', 'A8-08y12', 'A8-12y16', 'A8-16y20', 'A8-20y24', 
                          'A8-46Q66', 'A8-116Q150', 
                          'A13-norm',
                          'CMS7', 'CMS8', 
                          'CMS13-00y04','CMS13-04y08','CMS13-08y12','CMS13-12y16','CMS13-16y24',
                          #'CMS13_dQ_50to76',
                          'CMS13_dQ_106to170','CMS13_dQ_170to350','CMS13_dQ_350to1000',
                          'LHCb7', 'LHCb8', 'LHCb13_dy(2021)', 
                          'PHE200', 'STAR510', 
                          'E228-200', 'E228-300', 'E228-400', 
                          'E772',
                          'E605',
                          'D0run1-W','CDFrun1-W'
                          ]))

setDY=theData.CutData(cutFunc) 
#setDYfit=theData.CutData(cutFuncFORFIT) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')

#%%
import time
tt=[]
for i in range(25):
    t1=time.time()
    
    rSet.SetReplica(2*i+1)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
print("Avarage time of update:",numpy.mean(tt))

rSet.SetReplica(0)

#%%
DataProcessor.harpyInterface.PrintChi2Table(setDY,printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setSIDIS,printSysShift=False)

#%%
import pickle

with open(THIS_DIR+"/RESULT_1/"+"ART25","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2 = pickle.load(f)
    
X1=DataProcessor.harpyInterface.ComputeXSec(setDY)
X2=DataProcessor.harpyInterface.ComputeXSec(setSIDIS)

print((numpy.array(X1)-numpy.array(Y1))/numpy.array(X1))
print((numpy.array(X2)-numpy.array(Y2))/numpy.array(X2))

#%%
import pickle
ver="/RESULT_2/"

X1=DataProcessor.harpyInterface.ComputeXSec(setDY)
X2=DataProcessor.harpyInterface.ComputeXSec(setSIDIS)

# Saving the objects:
with open(THIS_DIR+ver+"ART25", 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([X1,X2], f)
    
#%%
import time
### Also save couple of replicas.
tt=[]
for i in [5,10,25,100]:
    t1=time.time()
    rSet.SetReplica(i)    
     
    X1=DataProcessor.harpyInterface.ComputeXSec(setDY)
    X2=DataProcessor.harpyInterface.ComputeXSec(setSIDIS)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
    
    with open(THIS_DIR+ver+"ART25_rep"+str(i), 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([X1,X2], f)
        
print("Avarage time of computation:",numpy.mean(tt))
#%%
ver="/RESULT_0/"
# Getting back the objects:
with open(THIS_DIR+ver+"snowART26","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2,Y3,Y4,Y5 = pickle.load(f)