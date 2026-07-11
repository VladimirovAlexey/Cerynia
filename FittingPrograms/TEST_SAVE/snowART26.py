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
import Cerynia
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

rSet=Cerynia.aTMDeReplicaSet.from_json(MODEL_DIR_ART25+"ART25_main.json")
rSet.set(0)

#%%

rSNOW=Cerynia.aTMDeReplicaSet.from_json(MODEL_DIR_SNOW26+"SnowART26.json")

rSNOW.set(0)
harpy.UpdateTables(1.0, 105.0)

#%%
def loadData(listOfNames, dirs):
    dataCollection = []
    for name in listOfNames:
        for d in dirs:
            path = DATAP_DIR + d + name + ".csv"
            if os.path.isfile(path):
                dataCollection.append(Cerynia.DataSet.from_csv(path))
                break
        else:
            raise FileNotFoundError(f"'{name}.csv' not found in any of {dirs}")
    return Cerynia.DataMultiSet(dataCollection)

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
def cutFuncD2(p):
    
    if p["Q_avg"]<numpy.sqrt(2.): return False, p
    
    if "E143" in p["id"] and p["proc_id"]==101:
        return False, p
    
    return True, p

##################Cut function
def cutFuncG2(p):
    if p["Q_avg"]<numpy.sqrt(2.):
        return False, p

    return True, p

##################Cut function
def cutFunc_tmdDY(p):    
    
    deltaTEST=0.3        
    delta=p["qT_avg"]/p["Q_avg"]        

    
    if(9<p["Q_avg"]<11):#UPSILON resonance-bin
        return False , p
    
    return delta<deltaTEST , p

def cutFunc_tmdSIDIS(p):     
    deltaTEST=0.35        
    if(p["pT_avg"]/p["z_avg"]<0.0):
        delta=0.0
        #delta=p["<pT>"]/p["<z>"]/p["<Q>"]        
    else:
        delta=p["pT_avg"]/p["z_avg"]/p["Q_avg"]        
        
    if "compass23" in p["id"]:
        deltaTEST=0.5
    if "compass08" in p["id"]:
        deltaTEST=0.47
    if "jlab" in p["id"]:
        deltaTEST=0.45
    if "JLab" in p["id"]:
        deltaTEST=0.45
   
    #### This is because star measures AN
    if "star" in p["id"]:
        p["thFactor"]=+p["thFactor"]
        
    return delta<deltaTEST and p["Q_avg"]>1.41, p


#%%
### Loading the SIDIS data set
theDataD2=loadData([
                "E143_d2",#"E154_d2",
                "E155-1999_d2","E155_d2",
                "HallA-2016_d2","HERMES_d2","SANE_d2",       
                #"RSS-2006_d2","RSS-2008_d2",    
                "RQCD_d2_ud",
                #"RQCD_d2_singlet","RQCD_d2_pn",
                #"GHMP26_d2","QCDSF_d2"
                ],
                ["DataLib/D2/"])

setD2=theDataD2.cut(cutFuncD2) 
setD2.prepare()

print(setD2)

#%%
### Loading the G2 data set
theDataG2=loadData([
    #"E142.n", 
    "E143.p", "E143.d",#"E143.n", 
    "E154.n",
    "E155-29.p","E155-32.p","E155-38.p",
    "E155-29.d","E155-32.d","E155-38.d",
    #"SMC.p",
    "HERMES",
    #"HallA-2004.n",
    "HallA-2016-4.He3","HallA-2016-5.He3"
    ],["DataLib/G2/"])

setG2=theDataG2.cut(cutFuncG2) 
setG2.prepare()

print(setG2)
#%%
### Loading the data set for Sivers
theDataS=loadData([
                    'compass08.sivers.pi+.dpt', 'compass08.sivers.pi-.dpt',
                    'compass08.sivers.k+.dpt', 'compass08.sivers.k-.dpt',
                    'compass08.sivers.pi+.dx', 'compass08.sivers.pi-.dx',
                    'compass08.sivers.k+.dx', 'compass08.sivers.k-.dx',
                    'compass08.sivers.pi+.dz', 'compass08.sivers.pi-.dz',
                    'compass08.sivers.k+.dz', 'compass08.sivers.k-.dz',
                    'compass16.sivers.h+.1<z<2.dpt','compass16.sivers.h-.1<z<2.dpt',
                    'compass16.sivers.h+.2<z.dpt' ,'compass16.sivers.h-.2<z.dpt',
                    'compass16.sivers.h+.1<z<2.dz','compass16.sivers.h-.1<z<2.dz',
                    'compass16.sivers.h+.2<z.dz' ,'compass16.sivers.h-.2<z.dz',
                    'compass16.sivers.h+.1<z<2.dx','compass16.sivers.h-.1<z<2.dx',
                    'compass16.sivers.h+.2<z.dx' ,'compass16.sivers.h-.2<z.dx',
                    'compass23.sivers.h+.dpt', 'compass23.sivers.h-.dpt',
                    'compass23.sivers.h+.dx', 'compass23.sivers.h-.dx',
                    'compass23.sivers.h+.dz', 'compass23.sivers.h-.dz',
                    'hermes3D.sivers.pi+','hermes3D.sivers.pi-',
                    'hermes3D.sivers.k+','hermes3D.sivers.k-',
                    'jlab.sivers.pi+','jlab.sivers.pi-','jlab.sivers.k+','jlab.sivers.k-'
                    ],["DataLib/Sivers/"])

theDataSdy=loadData([
                    'star26.sivers.W-.dy', 'star26.sivers.W+.dy',
                    'star23.sivers.Z'
                    ],["DataLib/Sivers/"])

setSivers=theDataS.cut(cutFunc_tmdSIDIS) 
setSivers.prepare()

setSiversDY=theDataSdy.cut(cutFunc_tmdDY) 
setSiversDY.prepare()

print(setSivers)
print(setSiversDY)



#%%
### Loading the WGT data set
theDataW=loadData([
                      'hermes3D.wgt.pi+','hermes3D.wgt.pi-',
                      'hermes3D.wgt.k+','hermes3D.wgt.k-',
                      'compass16.wgt.h+.2<z.dpt','compass16.wgt.h-.2<z.dpt',
                      'compass16.wgt.h+.2<z.dz','compass16.wgt.h-.2<z.dz',
                      'compass16.wgt.h+.2<z.dx','compass16.wgt.h-.2<z.dx',
                      'JLab6.wgt.pi+','JLab6.wgt.pi-'
                      ],["DataLib/wgt/"])

setALT=theDataW.cut(cutFunc_tmdSIDIS) 
setALT.prepare()

print(setALT)

#%%
weights_SiverDIS=Cerynia.harpyInterface.compute_weight(setSivers,method="central")
weights_SiverDY=Cerynia.harpyInterface.compute_weight(setSiversDY)
weights_wgt=Cerynia.harpyInterface.compute_weight(setALT,method="central")

#%%
Cerynia.harpyInterface.print_chi2_table(setD2,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setG2,sys_shift=False)

Cerynia.harpyInterface.print_chi2_table(setSivers,method="central",sys_shift=False,weights=weights_SiverDIS)
Cerynia.harpyInterface.print_chi2_table(setSiversDY,sys_shift=False,weights=weights_SiverDY)
Cerynia.harpyInterface.print_chi2_table(setALT,method="central",sys_shift=False,weights=weights_wgt)

#%%
import pickle

with open(THIS_DIR+"/RESULT_1/"+"snowART26","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2,Y3,Y4,Y5 = pickle.load(f)
    
X1=Cerynia.harpyInterface.xsec(setD2)
X2=Cerynia.harpyInterface.xsec(setG2)
X3=Cerynia.harpyInterface.xsec(setSivers,method="central",weights=weights_SiverDIS)
X4=Cerynia.harpyInterface.xsec(setSiversDY,weights=weights_SiverDY)
X5=Cerynia.harpyInterface.xsec(setALT,method="central",weights=weights_wgt)

print(numpy.array(X1)-numpy.array(Y1))
print(numpy.array(X2)-numpy.array(Y2))
print(numpy.array(X3)-numpy.array(Y3))
print(numpy.array(X4)-numpy.array(Y4))
print(numpy.array(X5)-numpy.array(Y5))

#%%
import pickle
ver="/RESULT_3/"

# Saving the objects:
with open(THIS_DIR+ver+"snowART26", 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([X1,X2,X3,X4,X5], f)
    
#%%
### Also save couple of replicas.
import time
tt=[]
for i in [5,10,25,100]:
    t1=time.time()
    rSNOW.set(i)
    rSet.set(rSNOW.get(i)['linkedReplica'][0])
    harpy.UpdateTables(1.0, 105.0) 
    
    weights_SiverDIS=Cerynia.harpyInterface.compute_weight(setSivers,method="central")
    weights_SiverDY=Cerynia.harpyInterface.compute_weight(setSiversDY)
    weights_wgt=Cerynia.harpyInterface.compute_weight(setALT,method="central")
     
    X1=Cerynia.harpyInterface.xsec(setD2)
    X2=Cerynia.harpyInterface.xsec(setG2)
    X3=Cerynia.harpyInterface.xsec(setSivers,method="central",weights=weights_SiverDIS)
    X4=Cerynia.harpyInterface.xsec(setSiversDY,weights=weights_SiverDY)
    X5=Cerynia.harpyInterface.xsec(setALT,method="central",weights=weights_wgt)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
    
    with open(THIS_DIR+ver+"snowART26_rep"+str(i), 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([X1,X2,X3,X4,X5], f)
        
print("Avarage time of computation:",numpy.mean(tt))