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


import sys
import numpy
sys.path.append(ROOT_DIR)
sys.path.append(HARPY_DIR)

###### RUN_CONFIGURATION to select what is what
#RUN_CONFIGURATION="LP"    ### LP usual
#RUN_CONFIGURATION="KPC0" ### KPC but with FUUT only
RUN_CONFIGURATION="KPC1" ### KPC includeing FUUL
#RUN_CONFIGURATION="KPC2" ### KPC including FUUL and boer-Collins function (TO BE MADE)

#%%
import Cerynia
import harpy


if RUN_CONFIGURATION=="LP" :
    path_to_constants=THIS_DIR+"/INI/ART25_main.atmde"
else:
    path_to_constants=THIS_DIR+"/INI/ART25_SIDIS_KPC.atmde"


harpy.initialize(path_to_constants)

inARRAY_TMDR=[1.5004, 0.073018, 0.038048, 0.0]
inARRAY_PDF=[0.521462, 0.000206, 0.402948, 7.0219, 1.0, 20.4051, 1.0, 0.000123, 1.1037, 0.660734, 0.0, 0.04]
inARRAY_FF=[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0]


harpy.setNPparameters_TMDR(inARRAY_TMDR)
harpy.setNPparameters_uTMDPDF(inARRAY_PDF)
harpy.setNPparameters_uTMDFF(inARRAY_FF)

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
### Loading the SIDIS data set
setSIDISc=loadData([
                'hermes3D.p.pi+','hermes3D.p.pi-',
                'hermes3D.d.pi+','hermes3D.d.pi-',
                'hermes3D.p.k+','hermes3D.p.k-',
                'hermes3D.d.k+','hermes3D.d.k-',
                'compass.d.h+','compass.d.h-'],
                ["/FittingPrograms/KPC26/DataCentered/"])

setSIDISc.prepare()

print(setSIDISc)

#%%
### Loading the DY data set
setDYc=loadData([
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
                ["/FittingPrograms/KPC26/DataCentered/"])

setDYc.prepare()

print(setDYc)
#%%

Cerynia.harpyInterface.print_chi2_table(setDYc,sys_shift=False,method="semiCentral")
Cerynia.harpyInterface.print_chi2_table(setSIDISc,sys_shift=False,method="semiCentral")

#%%
#######################################
# Minimisation
#######################################
import time

penalty_index=[-7,-6,-5,-4,-3]

def chi2(x):
    startT=time.time()
    #harpy.setNPparameters_uTMDFF([x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]])
    harpy.setNPparameters(x)
    print('np set =',["{:8.3f}".format(i) for i in x])        
    
    YY=Cerynia.harpyInterface.xsec(setSIDISc,method="semiCentral")
    ccSIDIS2,cc3=setSIDISc.chi2(YY)
        
    endTS=time.time()

    YY=Cerynia.harpyInterface.xsec(setDYc,method="semiCentral")
    ccDY2,cc3=setDYc.chi2(YY)
    
    
    #### This penalty term prevents low-energy DY to have extremely low normalization
    #penalty_array=numpy.array([max(0,abs(setDY.sets[i].DetermineAvarageSystematicShift(YY[setDY._i1[i]:setDY._i2[i]]))/setDY.sets[i].normErr[0]-1) for i in penalty_index])
    #penalty_term=sum(penalty_array**6)    
    
    ### This penalty term prevents SIDIS to be much lower than 1 (changes the slope of chi2 below 1)
    #pSIDIS=ccSIDIS2/setSIDIS.numberOfPoints
    #if pSIDIS<1.:
    #    penalty_term+=0.9*(1-pSIDIS)*setSIDIS.numberOfPoints    
    
    endT=time.time()
    
    print(':->    ',"{:.4f}".format((ccSIDIS2+ccDY2)/(setSIDISc.numberOfPoints+setDYc.numberOfPoints)),"    = ("
          " {:.2f}".format(ccSIDIS2/setSIDISc.numberOfPoints),
          " + ","{:.2f}".format(ccDY2/setDYc.numberOfPoints),          
          ")    time (DIS)=","{:.2f}".format(endTS-startT), 
          " time (DY)=","{:.2f}".format(endT-endTS)
          )
    
    return ccSIDIS2+ccDY2


#%%
#### Minimize SIDIS
from iminuit import Minuit

#---- PDFbias-like row (0.083931)
# initialValues=([1.5, 0.038969, 0.030641, 0.0, 
#                        0.51638, 0.002073, 0.478567, 0.373111, 
#                        2.407, 22.1996, 3.7876, 0.00128, 
#                        0.403343, 5e-05, 1.0, 1.0, 
#                        0.69769, 0.712969, -0.133895, -0.841651, 0.846846,
#                        0.774759, 1.5565, 1.1863, 0.692877, -0.569062, 
#                        0.0, 0.0])

initialValues=([1.5, 0.059952, 0.030591, 0.0, 
                       0.343842, 0.011714, 0.588161, 0.012262, 
                       7.682, 2.26267, 3.2397, 0.015181, 
                       0.052558, 1.2e-05, 1.0, 1.0, 
                       0.700419, 0.669345, -0.058065, -0.667958, 
                       0.863101, 0.7739, 1.396, 1.3471, 
                       0.654023, -0.490042, 0.0, 0.1])

initialErrors=(0.1,0.1,0.1,0.1,
                0.5,  1.0, 0.1,  1.0,
                0.5,  1.0, 0.1,  1.0,
                0.5,  1.0, 10.,  1.,
                0.5,0.5,0.5,0.5,
                0.5,0.5,0.5,0.5,
                0.5,0.5,0.5,0.5)
searchLimits=((1.0,2.5),(0.005,0.15) ,(0.0,.2), (-5.,5.),
              (0.00001,100.), (0.00001,100.),(0.00001,100.),(0.00001,100.),
              (0.00001,100.), (0.00001,100.),(0.00001,100.),(0.00001,100.),
              (0.00001,100.), (0.00001,100.),(0.0001,100.),(0.0001,100.),
              (0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100.),
              (0.0001,100.), (-100.,100.),(-100.,100.),(-100.,100),
              (-100.,100.),(-100.,100.),(-100.,100.),(0.1,2.))
              
# True= FIX
parametersToMinimize=(True, False,False,True,
                      False, False, False,False,
                      False, False, False, False,
                      False, False, True,True,
                      #True,True,True,True,
                      #True,True,True,True,
                      #True,True,True,True)
                      False, False, False,False,
                      False, False, False,False,
                      False, False, True,True)

#%%

m = Minuit(chi2, initialValues)

m.errors=initialErrors
m.limits=searchLimits
m.fixed=parametersToMinimize
m.errordef=1

print(m.params)
#%%
m.tol=0.0001*(setSIDISc.numberOfPoints+setDYc.numberOfPoints)*10000 ### the last 0.0001 is to compensate MINUIT def
m.strategy=1
m.migrad()

print(m.params)

chi2(list(m.values))
#%%
Cerynia.harpyInterface.print_chi2_table(setDYc,sys_shift=False,method="semiCentral")
Cerynia.harpyInterface.print_chi2_table(setSIDISc,sys_shift=False,method="semiCentral")

Cerynia.harpyInterface.print_chi2_table(setDYc,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setSIDISc,sys_shift=False)

print([round(x,1 if x >100 else 4 if x>1 else 6) for x in list(m.values)])

#%%
sys.exit()
