# HealthEconomics
Exit by strategy

David’s MATLAB scripts:

heOptimise.m
Runs the optimisation problem. Economic constraints hard coded in here. 

heRunCovid19.m
Use this to run a  simulation. Seasonality/seed etc set in here. Possibility of waning/drift and consecutive seasonal epidemics here also (turned off for now). 

heSimCovid19.m
Called by heRunCovid19 to run the core simulation. 

hePrepCovid19.m
Prepares as many parameters/objects as possible, including transmission matrix pre-lockdown. R0 defined in here and beta tuned accordingly. 

heMakeDs.m
Makes contact matrices for a given intervention stage. Input is population vector (by sector) and Xit (proportion of each sector open). The matrices A, B and C are hard coded in this script. 
