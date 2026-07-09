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

path_to_constants=THIS_DIR+"/INI/ART25_KT.atmde"


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
import time
tt=[]
for i in range(25):
    t1=time.time()
    
    rSet.SetReplica(2*i+1)
    
    t2=time.time()
    print("Update time =", t2-t1)
    tt.append(t2-t1)
print("Avarage time of update:",numpy.mean(tt))



#%%
import pickle

ver="/RESULT_1/"

xList=[0.0001,0.001,0.004,0.08,0.1,0.2,0.7]
zList=[0.06,0.1,0.2,0.7]
bList=[0.001,0.01,0.1,0.3,0.6,1.,1.5,2.,3.,5.]
kTList=[0.01,0.1,0.5,1.,2.,5.,10.,30.,100.]
QList=[-1.,2.,4.,10.,50.]


PDFlist_b=[]
for x in xList:
    for b in bList:
        for Q in QList:
            PDFlist_b.append(harpy.get_uTMDPDF(x, b, 1,mu=Q))
            
FFlist_b=[]
for z in zList:
    for b in bList:
        for Q in QList:
            FFlist_b.append(harpy.get_uTMDFF(z, b, 1,mu=Q))
            
PDFlist_k=[]
for x in xList:
    for kt in kTList:
        for Q in QList:
            PDFlist_k.append(harpy.get_uTMDPDF_kT(x, kt, 1,mu=Q))
            
FFlist_k=[]
for z in zList:
    for kt in kTList:
        for Q in QList:
            FFlist_k.append(harpy.get_uTMDFF_kT(z, b, 1,mu=Q))
#%%
with open(THIS_DIR+"/RESULT_0/"+"ART25_KT","rb") as f:  # Python 3: open(..., 'rb')
    PDFlist_b0,FFlist_b0,PDFlist_k0,FFlist_k0 = pickle.load(f)
    
print(numpy.array(PDFlist_b0)-numpy.array(PDFlist_b))
print(numpy.array(FFlist_b0)-numpy.array(FFlist_b))
print(numpy.array(PDFlist_k0)-numpy.array(PDFlist_k))
print(numpy.array(FFlist_k0)-numpy.array(FFlist_k))

#%%
# Saving the objects:
with open(THIS_DIR+ver+"ART25_KT", 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([PDFlist_b,FFlist_b,PDFlist_k,FFlist_k], f)