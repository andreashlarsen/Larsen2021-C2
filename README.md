# Larsen2021-C2

Files for generating simulations in the paper:     
Binding of Ca2+-Independent C2 Domains to Lipid Membranes: a Multi-Scale Molecular Dynamics Study
Andreas Haahr Larsen and Mark S. P. Sansom 

Published in Structure (open access):    
https://doi.org/10.1016/j.str.2021.05.011     

For questions/requests/collaborations, please contact Andreas Larsen: andreashlarsen@bioch.ox.ac.uk

## Overview of files:    

### Bash scripts   
LOGBOOK_C2.sh:        bash script to run coarse-grained simulations.   
UMBRELLA.sh:          bash script to run umbrella simulations and calculate PMFs.   
LOGBOOK_FEP:          bash script to run FEP calculations.    
analysis.sh:          bash script to analyze FEP calculations.    
CG2AT.sh:             bash script to convert from CG to AT with CG2AT.   
LOGBOOK_AT.sh:        bash script to run atomistic simulations.   
LOGBOOK_analyze_AT:   bash script to analyze atomistic simulations.    

### own python scripts
Dist_vs_time.py:      calculates and plots Dist vs. time (Fig. 3).       
Rzz_version2.py:      calculates and plots Rzz vs dist (Fig. 4).   
get_std_FEP.py:       calculates std dev for FEP.    
AT_analyze.py:        analyze and plot for AT sims (Fig. 6).    
Contacts.py:          calculate PIP2 contacts (Fig. 7).     
exponential_decay.py: fit exp docay to dist vs time (Fig. S1).     
pmf_vs_contacts.py:   plotting script for pmf vs contacts (Fig. S8).     
extract_frames:       extract frames for umbrella sampling, used in UMBRELLA.sh.      
Dist_PIP2.py:         
SHIP2_PMF_analyze:   
Dist_vs_time_PTEN:    calculate and plot dist vs time for larger-construct PTEN.      

### modified Martini scripts 
martinize_GROMACS_2018_plumed.py: edited line 1851 to have "/gromacs/top/" instead of "/top/"

### modified Martini topology files
martini_v2.0_lipids_all_201506.itp: added (optional) restaint on POP2, "POSRES_POP2"    

### plumed files
dist_ang_restraint_***.dat: Plumed files for each protein.   

### mdp files
##### note: some of the mdp files are copied and modified by LOGBOOK_C2.sh or UMBRELLA.sh or LOGBOOK_AT.sh before use.      
minimization.mdp:  CG minimization.   
equilibration.mdp: CG equilibration.    
production.mdp:    CG production.    
pull.mdp:          CG umbrella preparation.   
umbrella.mdp:      CG umbrella run.   
min.mdp:           AT minimization.   
nvt.mdp:           AT NVT equilibration.   
npt.mdp:           AT NPT equilibration.   
prod.mdp:          AT production.  
em_RFEP.mdp:       FEP equilibration.   
md_RFEP.mdp:       FEP production.     

### pdb files (structures) 
##### after modification by modeller.    
***.pdb.   
  
### Other files
LOGBOOK_Folders: text file with overview of the content in the different folders (referred to in the bash scripts).        
Overview.ods:    overview of C2 structures (used for initial selection of structures).   
