# Analysis of experiment e863
Code of the analysis of the experiment e863, a resonant elastic scattering experiment made at GANIL. The analysis is using [ROOT](https://root.cern/).

## To convert faster data :
```
cd faster-to-root/
```
_ Add the faster file in the **data/** folder 
_ Change the run number in the **faster2tree.C** code
```
#define DATAFILENAME "data/E863_run_0*run number*_0001.fast"
#define ROOTFILENAME "root/run_0*run number*.root"
```
_ Optional : check the type of data in the faster file 
```
faster_disfast -n *number of events* *name_of_the_file.fast*
```
_ To convert 
```
make clean
make
.\faster2tree.C
```
The output root file is saved in the **root/** folder.

## To lauch the post-analysis code : 
```
cd PostAnalysis/
root 'Runner.cxx("number_of_the_pipeline")'
```
Pipeline 1 : Applies energy and time calibration  
Pipeline 2 :  PID particle identification, DeltaE-TOF  
Pipeline 3 : Computes the energy un the center of mass

## Simulation : 
_ In **/sp_srim** stopping power tables from SRIM [J. F. Ziegler et al., SRIM Co., United States of America 6th (2013)] with S.P. in keV.micron  
_ loss_E_srim(Energy in MeV, distance in mm)  
_ Be careful while declaring the SRIM files to include the right number of lines in the file   
