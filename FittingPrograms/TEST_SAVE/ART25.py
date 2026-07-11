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
def loadData(listOfNames, dirs):
    dataCollection = []
    for name in listOfNames:
        for d in dirs:
            path = ROOT_DIR + d + name + ".csv"
            if os.path.isfile(path):
                dataCollection.append(Cerynia.DataSet.from_csv(path))
                break
        else:
            raise FileNotFoundError(f"'{name}.csv' not found in any of {dirs}")
    return Cerynia.DataMultiSet(dataCollection)

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
    
    if("1_W" in p["id"]):
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
theData=loadData([
                'hermes3D.p.pi+','hermes3D.p.pi-',
                'hermes3D.d.pi+','hermes3D.d.pi-',
                'hermes3D.p.k+','hermes3D.p.k-',
                'hermes3D.d.k+','hermes3D.d.k-',
                'compass.d.h+','compass.d.h-'],
                ["DataLib/SIDIS/"])

setSIDIS=theData.cut(cutFuncSIDIS) 
setSIDIS.prepare()

print(setSIDIS)

#%%
### Loading the DY data set
theData=loadData([
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
                ],
    
                ["DataLib/DY/","DataLib/DY_W/","DataLib/DY_angular/"])

setDY=theData.cut(cutFuncDY) 
setDY.prepare()

print(setDY)

#%%
rSet.set(0)
Cerynia.harpyInterface.print_chi2_table(setDY,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setSIDIS,sys_shift=False)

#%%
import pickle

with open(THIS_DIR+"/RESULT_2/"+"ART25","rb") as f:  # Python 3: open(..., 'rb')
    Y1,Y2 = pickle.load(f)
    
X1=Cerynia.harpyInterface.xsec(setDY)
X2=Cerynia.harpyInterface.xsec(setSIDIS)

print((numpy.array(X1)-numpy.array(Y1))/numpy.array(X1))
print((numpy.array(X2)-numpy.array(Y2))/numpy.array(X2))

#%%
import pickle
ver="/RESULT_3/"

X1=Cerynia.harpyInterface.xsec(setDY)
X2=Cerynia.harpyInterface.xsec(setSIDIS)

# Saving the objects:
with open(THIS_DIR+ver+"ART25", 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([X1,X2], f)

#%%
import time
tt=[]
for i in range(25):
    t1=time.time()
    
    rSet.set(2*i+1)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
print("Avarage time of update:",numpy.mean(tt))

rSet.set(0)
#%%
import time
### Also save couple of replicas.
tt=[]
for i in [5,10,25,100]:
    t1=time.time()
    rSet.set(i)
     
    X1=Cerynia.harpyInterface.xsec(setDY)
    X2=Cerynia.harpyInterface.xsec(setSIDIS)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
    
    with open(THIS_DIR+ver+"ART25_rep"+str(i), 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([X1,X2], f)
        
print("Avarage time of computation:",numpy.mean(tt))