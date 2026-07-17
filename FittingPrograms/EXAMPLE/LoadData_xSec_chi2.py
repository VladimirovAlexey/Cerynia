"""
Demonstrations of some elemental operations with Cerynia
In this example, it is shown
1) how to load data from the DataLib, and prepare it for computation
2) How to load artemide, and operate with replicas
3) how to compute cross-section and chi2 for given data
"""
#######################################
# The first part of any code is to load libraries
# In this example we need 
# 1) Cerynia (for operation with data)
# 2) harpy (to load artemide and operate with it)
#######################################

### Loading Cerynia is simple because we are already in its repository.
import Cerynia

### Loading harpy is a bit more cumbersome, because one should first specify the path to it.
### For future convenience I also save a couple more paths

import os
### this is just the current folder
THIS_DIR= os.path.dirname(__file__)


### This is the directory to the artemide
ATMDE_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..','..'))+"/artemide/"

#### all what is needed is to add harpy to PYTHONPATH
import sys
sys.path.append(ATMDE_DIR+"/harpy/")

#### Finally, I import harpy
import harpy

#%%
### To initialize artemide one needs to point to the INI-file, where the setup information is given
### NOTE: in this example I use ART25 model. It may not work (or work badly) with different setup
###       The model, and INI-file should be in agreement. In particular, the number of NP parameters of the model is
###       specified in the INI-file (and artemide will raise an error if one sends the wrong number). Also, different/specific
###       parts of artemide should be loaded to perform one-or-another computation. So, take care about this part.
###
### The following comand to harpy, starts artemide with a particular constant file
harpy.initialize(THIS_DIR+"/INI/ART25_main.atmde")

#%%
### At this moment artemide is initialized but it is not ready for computation.
### One should specify the values of NP parameters. It could be done in three ways
### 1) The principal way, is to call harpy.setNPparameters_???([list of floats]).
###     This will tell artemide the nonperturbative parameters for a particular input.
###     E.g. harpy.setNPparameters_TMDR([1.50, 0.073, 0.038, 0.0]), will set parameters of CS-kernel
###     NOTE: the meaning of parameters is specified in /src/Model/TMDR_model.f90 of artemide
###         In this particular example (ART25) it is:
###         DNP=NPparam(2)*bSS*bS+NPparam(3)*bSS*bS*Log(bS/NPparam(1))
###         So, the line harpy.setNPparameters_TMDR([1.50, 0.073, 0.038, 0.0]) sets
###         NPparam(1)=1.50, NPparam(2)=0.073, NPparam(3)=0.038.
###         The last 0.0 value is not used (and thus can be any number), but it is required because
###         the INI-file depends on a 4-parameter array.
###     Same can be done for other modules, such as _uTMDPDF, _uTMDFF, etc. 
###
### 2) The second way is to call harpy.setNPparameters([list of floats])
###    Here, [list]=[list for TMDR]+[list for uTMDPDF]+[list for uTMDFF]+....
###    In this case, artemide will split the list into parts and send each of them to necessary module
###    The order of elements is like the order of modules in INI-file
###
### 3) The last way, is to use an already saved list-of-NP-parameters (which is usually the result of a previous fit)
###    Usually, a list specifies a replica obtained in the fit, and thus we call it a replica.
###    A collection of replicas is saved in a json-file, and can be loaded into a Cerynia.aTMDeReplicaSet - object
###    In this case you should just specify the file, and the number of the replica. The rest is done automatically.
###
### Here I demonstrate the case 3).
### First, setup the aTMDeReplicaSet from the json (by ART25 replica-distribution, which is stored in ART25)

rSet=Cerynia.aTMDeReplicaSet.from_json(ATMDE_DIR+"Models/ART25/Replica-files/ART25_main.json")

### Then .set(n) will set the replica n.
### NOTE: n=0 is the mean replica
### NOTE: depending on setup it may take some time to set up a replica, because artemide may need to pre-compute tables
rSet.set(0)

###  At this moment artemide is ready for computation
### RECUP: To load artemide in python you need 3 steps
### 1) import harpy
### 2) harpy.initialize(path to INI-file)
### 3) Setup NP-parameters either, directly via harpy.setNPparameters or by saved-replicas rSet.set(n)

#%%
### SHORT DETOUR into aTMDeReplicaSet
### Using this class you can operate with replicas save/load jsons, set/read full and particular parts
### For example, .get(n) returns the dict of NP-parameters for replica n
print(rSet.get(0))
### .get(n,part="TMDR") returns only values for CS kernel
print(rSet.get(0,part="TMDR"))
### Or make a plot of distribution of parameters
import matplotlib.pyplot as plt
plt.hist([rSet.get(num)['TMDR']['params'][1] for num in range(1,rSet.numberOfReplicas+1)], bins=20)
plt.xlabel("c0")
plt.show()
### Or distributions of parameters of CS-kernel
plt.scatter([rSet.get(num)['TMDR']['params'][1] for num in range(1,rSet.numberOfReplicas+1)], 
           [rSet.get(num)['TMDR']['params'][2] for num in range(1,rSet.numberOfReplicas+1)],alpha=0.5)
#plt.colorbar(label="counts")
plt.xlabel("c0")
plt.ylabel("c1")
plt.show()

#%%
### Now, let us load the data
### The data can be created by yourself (see example CreateDataFile) or loaded from the library
### The data-library of Cerynia is located in /DataLib and it is a large collection of various experiments
### The synopsis of /DataLib/ is given in /DataLib/DataLib_Overview.xlsx
DATALIB_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..','..'))+"/DataLib/"
###
### The Cerynia has two classes related to Data-processing
### 1) DataSet represents the collection of points, that share same elements
### Such as common normalization, correlated uncertainty, covariance matrix, etc.
### 2) MultiDataSet, which is the collection of DataSets and is mainly needed to simplify computation
###  e.g. computation xSec for MultiDataSet is more efficient than the computation of xSec for each DataSet
### (due to proper parallelism)
###
### DataSets are saved in .csv files and can be loaded by Cerynia.DataSet.from_csv(path)
### Here, as an example, I load CDF1 measurement of Z-boson-production cross-section
CDF1_in=Cerynia.DataSet.from_csv(DATALIB_DIR+"DY/CDF1.csv")

### DataMultiSet can be made by adding together sets, or joining list of sets.
### Here, as an example, I load all E288 measurements (at 200, 300, and 400 GeVs)
E288_in=Cerynia.DataMultiSet([
    Cerynia.DataSet.from_csv(DATALIB_DIR+"DY/E288-200.csv"),
    Cerynia.DataSet.from_csv(DATALIB_DIR+"DY/E288-300.csv"),
    Cerynia.DataSet.from_csv(DATALIB_DIR+"DY/E288-400.csv")])

### You can inspect the sets by printing their information
print(CDF1_in)
print(E288_in)

#%%
### The DataSets contain all data points that are given by experimentalists.
### Naturally, we should apply some filtering of these points, e.g. kinematical cuts
### this is done by the command data.cut(cutFunc) , which returns a new DataSet/DataMultiSet
### which is obtained by applying cutFunc to each point.
###
### cutFunc is a function that receives the dict of a Point, and returns:
### the result of the test, and the point to be included
### NOTE: you can change the values of the point by cutFunc. It is very useful in analyses!
###
### Here is a simple function, which drops all points with 9<Q<11GeV (resonance region), and qT/Q>0.25 (out of TMD)
def cutFunc(p):
    if(9<p["Q_avg"]<11):#UPSILON resonance-bin
        return False , p   
    return p["qT_avg"]/p["Q_avg"]<0.25, p    
    
### It could be applied to DataSet or DataMultiSet
### NOTE: if the filter eliminates all points, it will return an empty dataset, which will be eliminated from DataMultiSet
CDF1=CDF1_in.cut(cutFunc)
E288=E288_in.cut(cutFunc)
### In this example, 27 points of CDF1 are dropped, and more than 200 from E288
print("In CDF1 it was: ", CDF1_in.numberOfPoints, " It became ",CDF1.numberOfPoints)
print("In E288 it was: ", E288_in.numberOfPoints, " It became ",E288.numberOfPoints)

#%%
### Finally, before operating with dataSets one should "prepare" them.
### This command computes such elements as covariance matrix, nuisance parameters, etc.
### NOTE: If a set is VERY large it may take some time, but usually it is very fast
CDF1.prepare()
E288.prepare()

### RECUP: To load and prepare data
### 1) import Cerynia
### 2) load DataSet/DataMultiSet from .csv file [by Cerynia.DataSet.from_csv(path)]
### 3) Apply filtering/cuts/modification by set.cut(OurFunction) [if needed]
### 4) set.prepare()  To compute the covariance matrix and other stuff.
###
### Now, we are all done!
### We have loaded and set up artemide, and we prepared the sets, you can start the analysis

#%%
### The main application is the computation of the cross-section for a given DataSet.
### It is done by Cerynia.harpyInterface.xsec(data)
### This command will send the data points to artemide, receive the result, and apply necessary transformations
### IMPORTANT NOTE: The artemide computes the average over the bin
###                 Particular definition depends on the process, see artemide-manual
###                 Meanwhile, each experiment presents something different.
###                 This factor is stored in p["thFactor"], which multiplies the result of computation of artemide
###                 and in this way gets agreement with experiment.
###                 Also, some experiments are normalized, .isNormalized=True. This is also taken into account in xsec()
###
### The result of computation is a list of floats for each point
X1=Cerynia.harpyInterface.xsec(CDF1)
X2=Cerynia.harpyInterface.xsec(E288)
### We can plot it. Since DataSet is a panda data-frame you can use its simplifications
plt.plot(CDF1.df["qT_avg"], X1, label="ART25", color="C0")
plt.errorbar(CDF1.df["qT_avg"], CDF1.df["xSec"], 
             yerr=(CDF1.df["uncorrErr_0"]), xerr=(CDF1.df["qT_max"]-CDF1.df["qT_min"])/2, 
             fmt="o", color="C1", label="CDF1", capsize=3)

plt.xlabel("qT")
plt.legend()
plt.show()

### There are several options of xsec()
### "method" specifies the way how computation is done (simplifications, speed-up, etc).
### It essentially changes the behavior, so read variants and their meaning in harpyInterface
### "weights" you can specify the weights to divide the results of computation,
### very useful for computation of asymmetries

#%%
### Another typical application of Cerynia is the computation of chi2
### This can be done by applying .chi2(X) to the values of cross-section.
### Such a command computes the chi2 of a given data set with given theoretical input.
### It also works for dataMultiSet, where it returns totalChi2, [list of chi2 for each exp.]
print(CDF1.chi2(X1))
print(E288.chi2(X2))

#%%
### As example, let make a histogram of chi2 for E288 fusing first 100 replicas
### This code takes some time to evaluate
ll=[]
for i in range(1,101):
    rSet.set(i) ### seting replica
    X2=Cerynia.harpyInterface.xsec(E288) ### compute cross-section
    ll.append(E288.chi2(X2)[0]/E288.numberOfPoints) ### make list of chi2/numPoints
plt.hist(ll, bins=15) ### Plot it
plt.xlabel("chi2 for E288")
plt.show()

#%%
### Finally, another useful function is printing the table of chi2
### It is done by the following command
### This command accepts "method" and "weights" from xsec, and also can print values of best-shift parameters
rSet.set(0)
Cerynia.harpyInterface.print_chi2_table(CDF1,sys_shift=True)
Cerynia.harpyInterface.print_chi2_table(E288,sys_shift=False)
