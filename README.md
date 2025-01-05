## PhD Thesis: Maximum Caliber approach to reweight dynamics of non-equilibrium steady states

[Link to Thesis](Thesis/thesis_final.pdf)

This is a collection of scripts used in my PhD thesis. Alot of the stuff here is purley experimental and not functionional without the proper setup, that I do not have anymore. 

I should have refactored/ cleaned up/ deleted / generalised after finishing the thesis. But I did not. So we have to live with what it is here. 

In case you want to recreate the contents of the thesis, I recommend following the detailed description in the thesis to implement your own code. You are also  welcome to contact me, so I can extract the important parts from the reopository. 

The workflow is like this: 
1. Simulate a model of your choice by  passing the paramaters: 


 (-arg (example value) description) 


```
/single_particle_model_sim/single_particle 

-lag (100 200 300 400) lagtime in units of muliple integreation steps dT. Multiple entries possible \
-lvl (6) number of energy levels \
-Ux (-U1 6 -U2 1 -U3 4 -U4 1 -U5 4 -U6 2 ) Energy level for each state \
-T (1) Temperature of heatbath\
-gamma (1) damping constant of heat bath\
-seed (0) seed of pRNG\
-pot (1) type of potential to be used  (1 or 2)\
-f (2) external force applied along x axis \
-wf (1/0) turn on/off writing trajectory details\
-fullT (1000)  simulation length\
-o $o (1000) output identifier\
-ms (30)  number of microstates\
-dT (0.00001) integration time step 
````

2. Perform markov check and MFPT/ fpt analysis:  

```
python MSM_analysis/check_markov.py \
-o (1000) input identififer \
-l (200) lagtime to analyse \
-rc (1) size of stable state in microstates\
-db (0) enforce detailed balance
```

3. Reweight two existing models into each other (for analysis) 

```
python reweighting/reweight_MSM_Sloc.py \
-i (1000) input identifier reweight from\
-o (2000) output identifier reweight to \
-l lagtime of MSM to reweight\
-rc size of stable stae in microstates 
```

