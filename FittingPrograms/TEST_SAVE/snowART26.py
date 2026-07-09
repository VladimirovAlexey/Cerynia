#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:54:42 2025

@author: alexey
"""

#######################################
# importing libraries
#######################################
import os
THIS_DIR= os.path.dirname(__file__)
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/"
DATAP_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))+"/"

SNOWFLAKE_DIR = ROOT_DIR+"artemide/harpy/"
MODEL_DIR_ART25 = ROOT_DIR+"artemide/Models/ART25/Replica-files/"
MODEL_DIR_SNOW26 = ROOT_DIR+"artemide/Models/SnowART26/Replicas/"

repFile=ROOT_DIR+"artemide/Models/SnowART26/Replicas/SnowART26"

import sys
import numpy
if('/data/arTeMiDe_Repository/artemide/harpy' in sys.path):
    sys.path.remove('/data/arTeMiDe_Repository/artemide/harpy')
sys.path.append(DATAP_DIR)
sys.path.append(SNOWFLAKE_DIR)

#%%
import DataProcessor.harpyInterface
import DataProcessor.snowInterface_N2
import DataProcessor.DataMultiSet
import harpy
#%%
#######################################
#Initialize snowflake
#######################################
path_to_INI=THIS_DIR+"/INI/snowflake_forRep.ini"
harpy.initialize_snowflake(path_to_INI)

NP_par=numpy.zeros(24)+0.2
harpy.setNPparameters_tw3(NP_par)
harpy.UpdateTables(1.0, 105.0)

#%%
#######################################
#Initialize artemide
#######################################

import DataProcessor.ArtemideReplicaSet

path_to_constants=THIS_DIR+"/INI/TMD+tw3.atmde"


harpy.initialize(path_to_constants)

rSet=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR_ART25+"ART25_main.rep")
    
rSet.SetReplica(0)

#%%
rSNOW=DataProcessor.ArtemideReplicaSet.ReadRepFile(MODEL_DIR_SNOW26+"SnowART26.rep")


rSNOW.SetReplica(0)
harpy.UpdateTables(1.0, 105.0)

#%%
### read the list of files and return the list of DataSets
def loadThisDataD2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/D2_moment/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection


### read the list of files and return the list of DataSets
def loadThisDataG2(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/G2/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   
        
    return dataCollection


### read the list of files and return the list of DataSets
def loadThisDataSivers(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/Sivers/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection


def loadThisDataWGT(listOfNames):    
    import DataProcessor.DataSet
    
    path_to_data=DATAP_DIR+"DataLib/wgt/"
    
    
    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)   

    return dataCollection
#%%
##################Cut function
def cutFunc(p):
    if p["type"]=="G2":
        if p["<Q>"]<numpy.sqrt(2.):
            return False, p
    
    return True, p

##################Cut function
def cutFuncD2(p):
    if p["type"]=="D2":
        if p["<Q>"]<numpy.sqrt(2.):
            return False, p
    
    if "E143" in p["id"] and p["process"]==101:
        return False, p
    
    return True, p

##################Cut function
def cutFunc_TMD(p):
    import copy
    
    if p["type"]=="DY":
        deltaTEST=0.3        
        delta=p["<qT>"]/p["<Q>"]        

        
        if(9<p["<Q>"]<11):#UPSILON resonance-bin
            return False , p
    
    if p["type"]=="SIDIS":   
        deltaTEST=0.35        
        if(p["<pT>"]/p["<z>"]<0.0):
            delta=0.0
            #delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
        else:
            delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
        
    if "compass23" in p["id"]:
        deltaTEST=0.5
    if "compass08" in p["id"]:
        deltaTEST=0.47
    if "jlab" in p["id"]:
        deltaTEST=0.45
    if "JLab" in p["id"]:
        deltaTEST=0.45
    
    if delta<deltaTEST:
        pNew=copy.deepcopy(p)    
        pNew["process"]=pNew["weightProcess"]
        if p["type"]=="SIDIS":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew,method="central")        
        elif p["type"]=="DY":
            normX=DataProcessor.harpyInterface.ComputeXSec(pNew)        
        else:
            print("Are you crazy?")
        p["thFactor"]=p["thFactor"]/normX        #### this minus is because of star
    
    #### This is because star measures AN
    if "star" in p["id"]:
        p["thFactor"]=+p["thFactor"]
    
    
        
#    return delta<0.5 and p.qT_avarage<80
    return delta<deltaTEST and p["<Q>"]>1.41, p


#%%
### Loading the D2 data set
theDataD2=DataProcessor.DataMultiSet.DataMultiSet("D2set",loadThisDataD2([
    "E143_d2",#"E154_d2",
    "E155-1999_d2","E155_d2",
    "HallA-2016_d2","HERMES_d2","SANE_d2",       
    #"RSS-2006_d2","RSS-2008_d2",    
    "RQCD_d2_ud",
    #"RQCD_d2_singlet","RQCD_d2_pn",
    #"GHMP26_d2","QCDSF_d2"
    ]))

setD2=theDataD2.CutData(cutFuncD2) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setD2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setD2.sets]), 'points.') 

#%%
### Loading the G2 data set
theDataG2=DataProcessor.DataMultiSet.DataMultiSet("G2set",loadThisDataG2([
    #"E142.n", 
    "E143.p", "E143.d",#"E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    "E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    "HERMES",
    #"HallA-2004.n",
    "HallA-2016-4.He3","HallA-2016-5.He3"
    ]))

setG2=theDataG2.CutData(cutFunc) 

#print('Loaded experiments are', [i.name for i in setDY.sets])

print('Loaded ', setG2.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setG2.sets]), 'points.') 

#%%
### Loading the data set for Sivers
theDataS=DataProcessor.DataMultiSet.DataMultiSet("SIDISset",loadThisDataSivers([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass08.sivers.pi+.dx', 'compass08.sivers.pi-.dx',
                    'compass08.sivers.k+.dx', 'compass08.sivers.k-.dx',
                    'compass08.sivers.pi+.dz', 'compass08.sivers.pi-.dz',
                    'compass08.sivers.k+.dz', 'compass08.sivers.k-.dz',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.z>2.dpt' ,'compass16.sivers.h-.z>2.dpt',
                    'compass16.sivers.h+.1<z<2.dz','compass16.sivers.h-.1<z<2.dz',
                    'compass16.sivers.h+.z>2.dz' ,'compass16.sivers.h-.z>2.dz',
                    'compass16.sivers.h+.1<z<2.dx','compass16.sivers.h-.1<z<2.dx',
                    'compass16.sivers.h+.z>2.dx' ,'compass16.sivers.h-.z>2.dx',
                    'compass23.sivers.h+.dpt', 'compass23.sivers.h-.dpt',
                    'compass23.sivers.h+.dx', 'compass23.sivers.h-.dx',
                    'compass23.sivers.h+.dz', 'compass23.sivers.h-.dz',
                    'hermes.sivers.pi+.3d','hermes.sivers.pi-.3d',
                    'hermes.sivers.k+.3d','hermes.sivers.k-.3d',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+','jlab.sivers.k-'
                    ]))

theDataSdy=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisDataSivers([
                    'star26.sivers.W-.dy', 'star26.sivers.W+.dy',
                    'star23.sivers.Z'
                    ]))

setSivers=theDataS.CutData(cutFunc_TMD) 
setSiversDY=theDataSdy.CutData(cutFunc_TMD) 


print('Loaded (SIDIS)', setSivers.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSivers.sets]), 'points.')
print('Loaded (DY)', setSiversDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setSiversDY.sets]), 'points.')
#print('Loaded SIDIS experiments are', [i.name for i in setSivers.sets])

#%%
### Loading the WGT data set
theDataW=DataProcessor.DataMultiSet.DataMultiSet("ALTset",loadThisDataWGT([
                      'hermes3D.ALT.pi+','hermes3D.ALT.pi-',
                      'hermes3D.ALT.k+','hermes3D.ALT.k-',
                      'compass16.ALT.h+.2<z.dpt','compass16.ALT.h-.2<z.dpt',
                      'compass16.ALT.h+.2<z.dz','compass16.ALT.h-.2<z.dz',
                      'compass16.ALT.h+.2<z.dx','compass16.ALT.h-.2<z.dx',
                      'JLab6.ALT.pi+','JLab6.ALT.pi-'
                      ]))

setALT=theDataW.CutData(cutFunc_TMD) 

print('Loaded ', setALT.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setALT.sets]), 'points.')

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setSiversDY,printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)
#%%
import time
tt=[]
for i in range(25):
    t1=time.time()
    
    rSNOW.SetReplica(2*i+1)
    rSet.SetReplica(rSNOW.GetReplica(2*i+1,part="ADD1")[0])
    harpy.UpdateTables(1.0, 105.0)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
print("Avarage time of update:",numpy.mean(tt))

rSNOW.SetReplica(0)
rSet.SetReplica(rSNOW.GetReplica(0,part="ADD1")[0])
harpy.UpdateTables(1.0, 105.0)

#%%
DataProcessor.snowInterface_N2.PrintChi2Table(setD2,printDecomposedChi2=False)
DataProcessor.snowInterface_N2.PrintChi2Table(setG2,printDecomposedChi2=False)

DataProcessor.harpyInterface.PrintChi2Table(setSivers,method="central",printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setSiversDY,printSysShift=False)
DataProcessor.harpyInterface.PrintChi2Table(setALT,method="central",printSysShift=False)

#%%
X1=DataProcessor.snowInterface_N2.ComputeXSec(setD2)
X2=DataProcessor.snowInterface_N2.ComputeXSec(setG2)
X3=DataProcessor.harpyInterface.ComputeXSec(setSivers,method="central")
X4=DataProcessor.harpyInterface.ComputeXSec(setSiversDY)
X5=DataProcessor.harpyInterface.ComputeXSec(setALT,method="central")

#%%
import pickle

with open(THIS_DIR+"/RESULT_1/"+"snowART26","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2,Y3,Y4,Y5 = pickle.load(f)
    
X1=DataProcessor.snowInterface_N2.ComputeXSec(setD2)
X2=DataProcessor.snowInterface_N2.ComputeXSec(setG2)
X3=DataProcessor.harpyInterface.ComputeXSec(setSivers,method="central")
X4=DataProcessor.harpyInterface.ComputeXSec(setSiversDY)
X5=DataProcessor.harpyInterface.ComputeXSec(setALT,method="central")

print(numpy.array(X1)-numpy.array(Y1))
print(numpy.array(X2)-numpy.array(Y2))
print(numpy.array(X3)-numpy.array(Y3))
print(numpy.array(X4)-numpy.array(Y4))
print(numpy.array(X5)-numpy.array(Y5))

#%%
import pickle
ver="/RESULT_2/"

# Saving the objects:
with open(THIS_DIR+ver+"snowART26", 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([X1,X2,X3,X4,X5], f)
    
#%%
### Also save couple of replicas.
import time
tt=[]
for i in [5,10,25,100]:
    t1=time.time()
    rSNOW.SetReplica(i)
    rSet.SetReplica(rSNOW.GetReplica(i,part="ADD1")[0])
    harpy.UpdateTables(1.0, 105.0) 
    
     
    X1=DataProcessor.snowInterface_N2.ComputeXSec(setD2)
    X2=DataProcessor.snowInterface_N2.ComputeXSec(setG2)
    X3=DataProcessor.harpyInterface.ComputeXSec(setSivers,method="central")
    X4=DataProcessor.harpyInterface.ComputeXSec(setSiversDY)
    X5=DataProcessor.harpyInterface.ComputeXSec(setALT,method="central")
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
    
    with open(THIS_DIR+ver+"snowART26_rep"+str(i), 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([X1,X2,X3,X4,X5], f)
        
print("Avarage time of computation:",numpy.mean(tt))
#%%
ver="/RESULT_0/"
# Getting back the objects:
with open(THIS_DIR+ver+"snowART26","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2,Y3,Y4,Y5 = pickle.load(f)

#%%
print(numpy.array(X1)-numpy.array(Y1))
print(numpy.array(X2)-numpy.array(Y2))
print(numpy.array(X3)-numpy.array(Y3))
print(numpy.array(X4)-numpy.array(Y4))
print(numpy.array(X5)-numpy.array(Y5))

#%%time
print(numpy.array(X1)-numpy.array(Y1))
