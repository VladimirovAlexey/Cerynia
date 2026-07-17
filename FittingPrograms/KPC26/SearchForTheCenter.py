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
rSet.set(0)
Cerynia.harpyInterface.print_chi2_table(setDY,sys_shift=False,method="semiCentral")
Cerynia.harpyInterface.print_chi2_table(setSIDIS,sys_shift=False,method="semiCentral")

#%%
import numpy as np

def searchBestPoint(ss,ls,tol):
    
    v=np.array([ss.df[r+"_avg"][0] for r in ls])
    num=int(1./tol)
    delta=np.array([(ss.df[r+"_max"][0]-ss.df[r+"_min"][0]) for r in ls])/2/num
    
    if(len(ls)==0):
        print("NOTHING TO DO")
        return ss
    if(len(ls)==1):
        vec=np.array([[+1],[-1]])
        kvec=[ls[0]+"_avg",ls[0]+"_avg"]
        ivec=[0,0]
    elif(len(ls)==2):
        vec=np.array([[+1,0],[-1,0],[0,+1],[0,-1]])
        kvec=[ls[0]+"_avg",ls[0]+"_avg",ls[1]+"_avg",ls[1]+"_avg"]
        ivec=[0,0,1,1]
    elif(len(ls)==3):
        vec=np.array([[+1,0,0],[-1,0,0],[0,+1,0],[0,-1,0],[0,0,+1],[0,0,-1]])
        kvec=[ls[0]+"_avg",ls[0]+"_avg",ls[1]+"_avg",ls[1]+"_avg",ls[2]+"_avg",ls[2]+"_avg"]
        ivec=[0,0,1,1,2,2]
    elif(len(ls)==4):
        vec=np.array([[+1,0,0,0],[-1,0,0,0],[0,+1,0,0],[0,-1,0,0],[0,0,+1,0],[0,0,-1,0],[0,0,0,+1],[0,0,0,-1]])
        kvec=[ls[0]+"_avg",ls[0]+"_avg",ls[1]+"_avg",ls[1]+"_avg",ls[2]+"_avg",ls[2]+"_avg",ls[3]+"_avg",ls[3]+"_avg"]
        ivec=[0,0,1,1,2,2,3,3]
    else:
        raise Exception("Only cases upt to 4 are made")
        
    f0=Cerynia.harpyInterface.xsec(ss)[0]
    fc=Cerynia.harpyInterface.xsec(ss,method="semiCentral")[0]-f0
    
    if(np.abs(fc/f0)<0.1*tolerance):
        print("IT IS ALREADY GOOD! (", f0, ":", fc+f0 ,")")
        return ss
    
    sign=+1 if fc>0 else -1
    fc=sign*fc
    
    #print(v,":",fc)
    for i in range(num):
        #### compute all values around given
        fMIN=fc
        jMIN=-1
        #print(jMIN,",",sign*fMIN+f0,",",f0," : ",fc)
        for j in range(len(vec)):
            ss.df.loc[0,kvec[j]]=(v+delta*vec[j])[ivec[j]]
            fCUR=sign*(Cerynia.harpyInterface.xsec(ss,method="semiCentral")[0]-f0)
            #print(j,",",sign*fCUR+f0,",",f0," : ",fCUR)
            if(fCUR<0):
                jMIN=j
                fMIN=fCUR
                break
            elif(fCUR<fMIN):
                jMIN=j
                fMIN=fCUR
        
        ###
        if(jMIN==-1): 
            print("THERE IS NO DESCEND! Return present point:", fCUR/f0)
            return ss
            
            
        #### now we know the direction of best case.
        #### May be it already the zero
        if(fMIN<0):
            v=v-delta*vec[jMIN]*fc/(fMIN-fc)
            for j in range(len(ls)): ss.df.loc[0,kvec[2*j]]=v[j]
            return ss
        else:
            #### go to next iteration
            v=v+delta*vec[jMIN]
            fc=fMIN
            for j in range(len(ls)): ss.df.loc[0,kvec[2*j]]=v[j]
            
            #print(v,":",fc)
            
        #if(i==num-1): raise Exception("AAA")
        if(i==num-1): print("DOES NOT FINISHED")
        
    return ss

#%%
tolerance=0.01
def applySearch(ss):
    if ss.processType=="DY":
        ll=[]
        if (ss.df["Q_max"][0]-ss.df["Q_min"][0])<=5.:
            ll.append("Q")
            ss.df.loc[0,"Q_avg"]=(ss.df.loc[0,"Q_min"]+ss.df.loc[0,"Q_max"])/2
            
        if (ss.df["qT_max"][0]-ss.df["qT_min"][0])<=4.:
            ll.append("qT")
            ss.df.loc[0,"qT_avg"]=(ss.df.loc[0,"qT_min"]+ss.df.loc[0,"qT_max"])/2
            
        if (ss.df["y_max"][0]-ss.df["y_min"][0])<=2.:
            ll.append("y")           
            ss.df.loc[0,"y_avg"]=(ss.df.loc[0,"y_min"]+ss.df.loc[0,"y_max"])/2
        
        ss.isNormalized=False
        #### first approximate search
        ss= searchBestPoint(ss,ll,0.25)
        #### than accurate
        return searchBestPoint(ss,ll,tolerance)
    
    if ss.processType=="SIDIS":
        ll=[]
        if (ss.df["Q_max"][0]-ss.df["Q_min"][0])<5.:
            ll.append("Q")
            ss.df.loc[0,"Q_avg"]=(ss.df.loc[0,"Q_min"]+ss.df.loc[0,"Q_max"])/2
            
        if (ss.df["pT_max"][0]-ss.df["pT_min"][0])<4.:
            ll.append("pT")
            ss.df.loc[0,"pT_avg"]=(ss.df.loc[0,"pT_min"]+ss.df.loc[0,"pT_max"])/2
            
        ll.append("x")
        ll.append("z")
        ss.df.loc[0,"x_avg"]=(ss.df.loc[0,"x_min"]+ss.df.loc[0,"x_max"])/2
        ss.df.loc[0,"z_avg"]=(ss.df.loc[0,"z_min"]+ss.df.loc[0,"z_max"])/2
        ss= searchBestPoint(ss,ll,0.25)
        return searchBestPoint(ss,ll,tolerance)
    return ss
#%%
for s in setDY:
    new_ds=Cerynia.DataSet.empty(s.processType, name=s.name, normErr=s.normErr,
                          isNormalized=s.isNormalized, normalizationMethod=s.normalizationMethod)
    for i in range(s.numberOfPoints):
        print(new_ds.name,":",i,"/",s.numberOfPoints)
        sp=applySearch(s.point(i))
        new_ds.add_point(sp.df.iloc[0].to_dict())    
    new_ds.save_csv(THIS_DIR+"/DataCentered/" + new_ds.name + ".csv")

#%%
for s in setSIDIS:
    new_ds=Cerynia.DataSet.empty(s.processType, name=s.name, normErr=s.normErr,
                          isNormalized=s.isNormalized, normalizationMethod=s.normalizationMethod)
    for i in range(s.numberOfPoints):
        print(new_ds.name,":",i,"/",s.numberOfPoints)
        sp=applySearch(s.point(i))
        new_ds.add_point(sp.df.iloc[0].to_dict())        
    new_ds.save_csv(THIS_DIR+"/DataCentered/" + new_ds.name + ".csv")

exit()
#%%
for i in range(5):
    SS=setDY.point(350+i)
    print(SS)
    W=applySearch(SS)
    #SS.df.loc[0,"qT_avg"]=SS.df.loc[0,"qT_max"]
    X1=Cerynia.harpyInterface.xsec(W)
    X2=Cerynia.harpyInterface.xsec(W,method="semiCentral")
    print("{",W.df['qT_avg'][0],",",W.df['y_avg'][0],",",X1[0],",",X2[0],"},")

#%%
SS=setDY.point(120)
#W=applySearch(SS)
for i in range(10):
    SS.df.loc[0,"qT_avg"]=SS.df.loc[0,"qT_min"]+(SS.df.loc[0,"qT_max"]-SS.df.loc[0,"qT_min"])*i/10
    W=SS
    X1=Cerynia.harpyInterface.xsec(W)
    X2=Cerynia.harpyInterface.xsec(W,method="semiCentral")
    print("{",W.df['qT_avg'][0],",",W.df['y_avg'][0],",",X1[0],",",X2[0],"},")

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
Cerynia.harpyInterface.print_chi2_table(setDY,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setSIDIS,sys_shift=False)

Cerynia.harpyInterface.print_chi2_table(setDYc,sys_shift=False)
Cerynia.harpyInterface.print_chi2_table(setSIDISc,sys_shift=False)
#%%
Cerynia.harpyInterface.print_chi2_table(setDYc,sys_shift=False,method="semiCentral")
Cerynia.harpyInterface.print_chi2_table(setSIDISc,sys_shift=False,method="semiCentral")