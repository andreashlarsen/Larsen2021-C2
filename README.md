# Larsen2021-C2

Files for generating simulations in the paper:     
Binding of Ca2+-Independent C2 Domains to Lipid Membranes: a Multi-Scale Molecular Dynamics Study
Andreas Haahr Larsen and Mark S. P. Sansom 

BioRxiv:      
https://www.biorxiv.org/content/10.1101/2020.10.30.361964v1.   

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
AT_analyze.py:        analyze and plot for AT sims (Fig. 6). 
Contacts.py:          calculate PIP2 contacts (Fig. 7).  
exponential_decay.py: fit exp docay to dist vs time (Fig. S1)
pmf_vs_contacts.py:   plotting script for pmf vs contacts (Fig. S8).
extract_frames:       extract frames for umbrella sampling, used in UMBRELLA.sh.   
Dist_PIP2.py:         
SHIP2_PMF_analyze:   
Dist_vs_time_PTEN:    calculate and plot dist vs time for larger-construct PTEN.     

### modified Martini scripts 
martinize_GROMACS_2018_plumed.py: edited line 1851 to have "/gromacs/top/" instead of "/top/"

### plumed files
dist_ang_restraint_<XXX>.dat: Plumed files for each protein.   

### mdp files
#### note: some of the mdp files are copied and modified by LOGBOOK_C2.sh or UMBRELLA.sh or LOGBOOK_AT.sh before use.      
minimization.mdp:  CG minimization.   
equilibration.mdp: CG equilibration.    
production.mdp:    CG production.    
pull.mdp:          CG umbrella preparation.   
umbrella.mdp:      CG umbrella run.   
min.mdp:           AT minimization.   
nvt.mdp:           AT NVT equilibration.   
npt.mdp:           AT NPT equilibration.   
prod.mdp:          AT production.  

### pdb files (structures) after modification by modeller.  
<XXX>.pdb.   
  
### Other files
LOGBOOK_Folders: text file with overview of the content in the different folders (referred to in the bash scripts).   
Overview.ods:    overview of C2 structures (used for initial selection of structures).   
