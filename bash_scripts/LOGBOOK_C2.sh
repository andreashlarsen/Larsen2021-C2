#!/bin/sh
#

# --------------------------------------------------------------
# 
# This bash script runs simulations with C2 domains from different proteins and the following membranes
#
# C2 + POPC:POPS:PIP2:PIP3 (80:13:5:2) Bilayer
#
# C2 + POPC:POPS:PIP2(80:15:5) Bilayer
#
# C2 + POPC:POPS (80:20) Bilayer
#
# C2 + POPC Bilayer
# 
# CG MARTINI 2.1
# 
# Andreas Haahr Larsen
#
# November 2019
# 
# --------------------------------------------------------------

## ensure modules and other settings from batchrc file
source ~/.bashrc

## set parameters

###################### expect user input from here ##################################

# select protein (8=PI3K, 9=FVa,11=RIM2, 12=KIBRA, 13=MFGE8, 14=SHIP2, 15=SMURF2, 16=PTEN, 17=FVIII, 18=RIM2_t, 19=DLL1, 20=DLL4, 21=DLL4_t, 22=KIBRA_t, 23=PTEN_t, 24=SHIP2_t, 24=SHIP2_tt, 28=SHIP2_FL, 29=PTEN_FL)
protein_folder=15

# set lipids to loop over (1=PC:PS:PIP2:PIP3, 2=PC:PS:PIP2, 3=PC:PS, 4=PC, 5=CHOL:PC:PS:PIP2(25:55:15:5), 6=PC:PS (60:40), 7=PC:PS:PIP1, 8=PC:PS:PIP3, 9=CHOL:PC:PS(25:55:20), 10=CHOL:PC(25:75))
lip_first=0
lip_last=2

# set number of repetitions to run for each lip
rep_first=19
rep_last=24

# activate/deactivate overall modules of script
PREPARE=0
AXON=0
MDRUN=0
DOWNLOAD=0
ANALYSIS=1

# submodules of ANALYSIS (only relevant if ANALYSIS=1)
ANALYSIS_CENT=0
ANALYSIS_DIST=0
ANALYSIS_DIST_COM=0
ANALYSIS_DENS=0
ANALYSIS_RZZ=1
ANALYSIS_RMSD=0
ANALYSIS_SASA=0
ANALYSIS_CONTACT=0
ANALYSIS_MSD=0

# parallelisation with mdrun
pinoffset=0
noCPUs=2
noGPUs=1

# nb=auto to use gpu, nb=cpu to avoid use of gpu (if used by other jobs)
nb=auto

########################## to here ###################################################

# set parameters depending on protein group
NEUTRAL_PROTEIN=0 # by default assume that protein has some net charge
if [ $protein_folder -eq 8 ]
then
  protein_name=PI3K
  pdb_name=6btz
  lip_for_ref=2
  key_frame_rep=1
  key_frame=1420
elif [ $protein_folder -eq 9 ]
then
  protein_name=FV 
  pdb_name=1czv
  lip_for_ref=2
  key_frame_rep=18
  key_frame=4379
elif [ $protein_folder -eq 11 ]
then
  protein_name=RIM2
  pdb_name=2bwq
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 12 ]
then
  protein_name=KIBRA
  pdb_name=6fjd
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 13 ]
then
  protein_name=MFGE8

  key_frame_rep=21
  key_frame=2039
elif [ $protein_folder -eq 14 ]
then
  protein_name=SHIP2
  pdb_name=5okm
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 15 ]
then
  protein_name=SMURF2
  pdb_name=2jqz # 2jqz_ss is with secondary structure preserved
  lip_for_ref=2
  key_frame_rep=24
  key_frame=2299
elif [ $protein_folder -eq 16 ]
then
  protein_name=PTEN 
  pdb_name=1d5r
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 17 ]
then
  protein_name=FVIII 
  pdb_name=3hny
  lip_for_ref=2
  key_frame_rep=6
  key_frame=5594
elif [ $protein_folder -eq 18 ]
then
  protein_name=RIM2
  pdb_name=2bwq_t
  lip_for_ref=2
  key_frame_rep=11
  key_frame=4849
elif [ $protein_folder -eq 19 ]
then
  protein_name=DLL1
  pdb_name=4xbm
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 20 ]
then
  protein_name=DLL4
  pdb_name=5mvx
  lip_for_ref=2
  key_frame_rep=1
  key_frame=1
  NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 21 ]
then
  protein_name=DLL4
  pdb_name=5mvx_t
  lip_for_ref=2
  key_frame_rep=15
  key_frame=5483
  NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 22 ]
then
  protein_name=KIBRA 
  pdb_name=6fjd_t
  lip_for_ref=2
  key_frame_rep=11
  key_frame=4665
elif [ $protein_folder -eq 23 ]
then
  protein_name=PTEN 
  pdb_name=1d5r_t
  lip_for_ref=2
  key_frame_rep=1
  key_frame=2431
elif [ $protein_folder -eq 24 ]
then
  protein_name=SHIP2
  pdb_name=5okm_t
  lip_for_ref=2
  key_frame_rep=0
  key_frame=0
elif [ $protein_folder -eq 25 ]
then
  protein_name=SHIP2
  pdb_name=5okm_tt
  lip_for_ref=2
  key_frame_rep=23
  key_frame=4932
elif [ $protein_folder -eq 26 ]
then
  protein_name=DLL1
  pdb_name=4xbm_t
  lip_for_ref=2
  key_frame_rep=2
  key_frame=2246
elif [ $protein_folder -eq 27 ]
then
  protein_name=PI3K
  pdb_name=6bu0
  lip_for_ref=2
  key_frame_rep=7
  key_frame=5108
elif [ $protein_folder -eq 28 ]
then
  protein_name=SHIP2_FL
  pdb_name=5okm_FL
  lip_for_ref=2
  key_frame_rep=13
  key_frame=5611
elif [ $protein_folder -eq 29 ]
then
  protein_name=PTEN_FL
  pdb_name=1d5r_FL
  lip_for_ref=2
  key_frame_rep=3
  key_frame=5000
elif [ $protein_folder -eq 99 ] # for testing purpose
then
  protein_name=TEST # protein name 
  pdb_name=test # # pdb name 
  lip_for_ref=2 # lipid combi used to generate reference structure (for Rzz analysis). 1 = PC:PS:PIP2:PIP3, 2 = PC:PS:PIP2, 3 = PC:PS, 4 = PC
  key_frame_rep=0 # repetition of frame used to generate reference structure (for Rzz analysis)
  key_frame=0 # number of frame used to generate reference structure (for Rzz analysis)
  NEUTRAL_PROTEIN=1
fi

# define paths to scripts 
py2=/sansom/s157/bioc1642/anaconda2/bin/python2.7
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
martinize=/sansom/s157/bioc1642/Desktop/Scripts/martinize_GROMACS_2018_plumed.py # edited line 1851 to have "/gromacs/top/" instead of "/top/"
insane=/sansom/s157/bioc1642/Desktop/Scripts/insane.py

# define paths to mdp files
min=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/minimization.mdp
eq=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/equilibration.mdp
prod=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp
 
# define path to pdb file
protein=/sansom/s157/bioc1642/Desktop/prj_C2/Structures/${pdb_name}_clean.pdb # truncated

# define path to topology files
ffdir=/sansom/s157/bioc1642/Desktop/Scripts/martini2

# plumed files
plumed_file=/sansom/s157/bioc1642/Desktop/Scripts/plumed/dist_ang_restraint_${protein_name}.dat

# simulation time in ps to start collecting density maps
time=1996000

## create and go to new protein working directory
mkdir -p $protein_folder
cd $protein_folder

## loop over lipid compositions
for j in $(seq $lip_first $lip_last)
do

  # folder name
  if [ $j -eq 1 ]
  then
    Folderprefix=PCPSPIP2PIP3
    short_lip_name=P3
  elif [ $j -eq 2 ]
  then
    Folderprefix=PCPSPIP2
    short_lip_name=PI
  elif [ $j -eq 3 ]
  then
    Folderprefix=PCPS
    short_lip_name=PS
  elif [ $j -eq 4 ]
  then
    Folderprefix=PC
    short_lip_name=PC
  elif [ $j -eq 5 ]
  then
    Folderprefix=CHPCPSPIP2
    short_lip_name=CH
  elif [ $j -eq 6 ]
  then
    Folderprefix=PCPS40
    short_lip_name=PS40
  elif [ $j -eq 7 ]
  then  
    Folderprefix=PCPSPIP1
    short_lip_name=P1
  elif [ $j -eq 8 ]
  then  
    Folderprefix=PCPSPIP3
    short_lip_name=P3  
  elif [ $j -eq 9 ]
  then
    Folderprefix=CHPCPS
    short_lip_name=CHPS 
  elif [ $j -eq 10 ]
  then
    Folderprefix=CHPC
    short_lip_name=CHPC
  fi

  ## loop over replicas
  for i in $(seq $rep_first $rep_last)
  do
    
    ############ GENERATE FOLDERS ETC #######################################################
    
    ## create (if not already existing) and go to working directory
    folder=${Folderprefix}_$i
    mkdir -p $folder
    cd $folder
   
    ## print directory
    echo ----------------------------
    echo $folder
    echo ----------------------------- 
    ############ PREPARE TO RUN SIMULATION ###################################################
    if [ $PREPARE -eq 1 ]
    then
      ## rotate protein 
      gmx insert-molecules -ci $protein -o ${protein_name}_rot.pdb -nmol 1 -box 6 6 6 -rot xyz -quiet

      ## center protein
      gmx editconf -f ${protein_name}_rot.pdb -o ${protein_name}_rot_cent.pdb -c -quiet

      ## CG protein  with martinize 2.2
      $py3 $martinize -f ${protein_name}_rot_cent.pdb -x ${protein_name}_CG.pdb -o topol.top -v -ff martini22 -elastic -dssp dssp

      ## Build bilayer + protein with insane with 10% antifreeze
      box_xy=7
      box_z=18	
      if [ $protein_folder -eq 28 ] || [ $protein_folder -eq 29 ]
      then
        box_xy=12
        box_z=24
	echo "I am here"
      fi
      dist=4.4 # four times cutoff length (see mpd file)
      if [ $j -eq 1 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:13 -l POP2:5 -l POP3:2 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 2 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:15 -l POP2:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist 
      elif [ $j -eq 3 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:20 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 4 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 5 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l CHOL:25 -l POPC:55 -l POPS:15 -l POP2:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 6 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:60 -l POPS:40 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 7 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:15 -l POP1:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 8 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:80 -l POPS:15 -l POP3:5 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      elif [ $j -eq 9 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l CHOL:25 -l POPC:55 -l POPS:20 -sol W:9 -sol WF:1 -salt 0 -dm $dist      
      elif [ $j -eq 10 ]
      then
        $py2 $insane -f ${protein_name}_CG.pdb -o ${protein_name}_Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l CHOL:25 -l POPC:75 -sol W:9 -sol WF:1 -salt 0 -dm $dist
      fi

      ## change typology file for system
      sed -i -e 's/#include "martini.itp"//g' topol.top
      cat << EOF > topol.add
#include "$ffdir/martini_v2.2.itp"
#include "$ffdir/martini_v2.0_ions.itp"
#include "$ffdir/martini_v2.0_lipids_all_201506.itp"
#include "Protein.itp"

#ifdef POSRES
#include "posre.itp"
#endif
EOF
      cat topol.add topol.top > tmp
      rm topol.add
      mv tmp topol.top

      ## minimization
      gmx grompp -f $min -c ${protein_name}_Bilayer.gro -p topol.top -o min.tpr -quiet
      gmx mdrun -deffnm min -quiet -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb

      ## make index file
      if [ $j -eq 1 ]
      then
        cat << EOF > index.input
17 | 18 | 19 | 25 | 26 | 27
name 28 SOL
13 | 14 | 15 | 16 | 21 | 22 | 23 | 24
name 29 LIP
a BB
14 | 22
name 31 PS
15 | 23
name 32 PIP2
16 | 24
name 33 PIP3
q
EOF
      elif [ $j -eq 2 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 20
name 29 PS
15 | 22
name 30 PIP2
q
EOF
      elif [ $j -eq 3 ]
      then 
        cat << EOF > index.input
15 | 16 | 17 | 21 | 22 | 23
name 24 SOL
13 | 14 | 19 | 20
name 25 LIP
a BB
14 | 21
name 27 PS
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 1 ]
      then
        cat << EOF > index.input
14 | 15
name 16 SOL
13
name 17 LIP
a BB
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 0 ]
      then
        cat << EOF > index.input
14 | 15 | 16 | 19 | 20 | 21
name 22 SOL
13
name 23 LIP
a BB
q
EOF
      elif [ $j -eq 5 ]
      then
        cat << EOF > index.input
17 | 18 | 19 | 25 | 26 | 27
name 28 SOL
13 | 14 | 15 | 16 | 21 | 22 | 23 | 24
name 29 LIP
a BB
15 | 23
name 31 PS
16 | 24
name 32 PIP2
q
EOF
      elif [ $j -eq 6 ]
      then
        cat << EOF > index.input
15 | 16 | 17 | 21 | 22 | 23
name 24 SOL
13 | 14 | 19 | 20
name 25 LIP
a BB
14 | 20
name 27 PS
q
EOF
      elif [ $j -eq 7 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 21
name 29 PS
15 | 22
name 30 PIP1
q
EOF
      elif [ $j -eq 8 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
14 | 21
name 29 PS
15 | 22
name 30 PIP3
q
EOF
      elif [ $j -eq 9 ]
      then
        cat << EOF > index.input
16 | 17 | 18 | 23 | 24 | 25
name 26 SOL
13 | 14 | 15 | 20 | 21 | 22
name 27 LIP
a BB
15 | 22
name 29 PS
q
EOF
      elif [ $j -eq 10 ]
      then
        cat << EOF > index.input
15 | 16 | 17 | 21 | 22 | 23
name 24 SOL
13 | 14 | 19 | 20
name 25 LIP
a BB
q
EOF
      fi
      gmx make_ndx -f ${protein_name}_Bilayer.gro -quiet < index.input
      rm index.input

      ## see overview of index file 
      # gmx make_ndx -f ${protein_name}_Bilayer.gro -n index.ndx -quiet
  
      ## generate position restraint for protein
      echo 1 | gmx genrestr -f ${protein_name}_CG.pdb -fc 1000 1000 1000 -o posre.itp -quiet

      ## equilibrate, protein restrained, wall restraint
      if [ $j -le 2 ] 
      then
        cp $eq eq_dt.mdp
        sed -i -e 's/dt                       = 0.03/dt                       = 0.02/g' eq_dt.mdp # lower dt from 0.03 to 0.02      
        sed -i -e 's/nsteps                   = 300000/nsteps                   = 500000/g' eq_dt.mdp
        sed -i -e 's/define                   = -DSTRONG_POSRES   ; Prevent protein from moving too much/define                   = -DSTRONG_POSRES -DFLEXIBLE   ; Prevent protein from moving too much and change some constraints in PIP to bonds/g' eq_dt.mdp # change some constrants in pip to bonds
        gmx grompp -f eq_dt.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -quiet -n index.ndx
      else
        gmx grompp -f $eq -c min.gro -r min.gro -p topol.top -o eq.tpr -quiet -n index.ndx
      fi
      gmx mdrun -deffnm eq -v -quiet -plumed $plumed_file -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb
      mv COLVAR COLVAR_equilibration

      ## prepare production run
      cp $prod production_short.mdp
      sed -i -e 's/nsteps                   = 86000000/nsteps                   = 57333333/g' production_short.mdp # lower production time from 3 to 2 us
      if [ $protein_folder -eq 28 ]
      then
        sed -i -e 's/tau_p                    = 14.0/tau_p                    = 24.0/g' production_short.mdp
      fi
      gmx grompp -f production_short.mdp -c eq.gro -r eq.gro -p topol.top -o md.tpr -quiet -n index.ndx
    
    # end PREPARE if statement
    fi
    
    ############ RUN SIMULATION ##############################################################
    if [ $MDRUN -eq 1 ]
    then
    	gmx mdrun -deffnm md -v -quiet -plumed $plumed_file -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb
    fi
    
    ############ COPY TO AXON #################################################################
    if [ $AXON -eq 1 ]
    then
	cat << EOF > Flow_Axon.sh
#!/bin/bash 
## Set the job name.
#SBATCH --job-name=${protein_folder}_${short_lip_name}_$i

## Set the number of nodes and cores. This shouldn't need changing.
## The maximum number of nodes which can be used on axon is 1.
#SBATCH --nodes=1

## Set the number of tasks on each node, this is usually the number of program executions done.
## For example, if you are running an mpi run, this would reflect the number passed to the -np flag.
#SBATCH --ntasks-per-node=1

## Set the number of cores per tasks. This usually reflects the number of threads (often openmp) that
## are being assigned per task. Benchmarking has shown that setting tasks to 1 and cpus-per-task to
## 16 for atomistic simulations and 8 for CG is a good starting point for getting maximum efficiency.
#SBATCH --cpus-per-task=6

## Set the number of GPUs to be used. In most cases this will be set to 1.
#SBATCH --gres=gpu:1

## IMPORTANT: set GPU binding, otherwise your jobs will clash with other users'
#SBATCH --gres-flags=enforce-binding

## Select the queues you will be running on (sansom: gpu-sansom,gpu-sansom2 biggin: gpu-biggin,gpu-biggin2) 
##SBATCH -p gpu-sm-short
##SBATCH -p gpu-sm-urgent
#SBATCH -p gpu-sansom, gpu-sansom2

## Select the max amount of time this job will run (48h for gpu-sansom, 3h for shor)
##SBATCH --time=3:00:00
#SBATCH --time=48:00:00

source /etc/profile.d/modules.sh
module purge
module load apps/gromacs/2018.6-plumed_2.4.4-GPU-KEPLER

## Note if running an energy minimisation, CG or using energy groups, you need to add the -ntmpi 1 for gromacs 2019 and above
## production, continuation (500000 = 1 ns, 50000000 = 100 ns, 100000000 = 200 ns)
gmx mdrun -deffnm md -v -quiet -plumed plumed.dat -ntomp \${SLURM_CPUS_PER_TASK} -ntmpi 1 -rdd 1.6
EOF

        axon_folder=${protein_name}_${short_lip_name}_$i
        mkdir -p $axon_folder
        cp md.tpr $axon_folder
        cp $plumed_file $axon_folder/plumed.dat
	mv Flow_Axon.sh $axon_folder
        scp -r $axon_folder axon:/home/bioc1642/
    fi


    ############ DOWNLOAD FROM AXON ###########################################################
   
    if [ $DOWNLOAD -eq 1 ]
    then
      axon_folder=${protein_name}_${short_lip_name}_$i
      scp axon:/home/bioc1642/${axon_folder}/md.* .
    fi
    ############ ANALYSIS #####################################################################
    if [ $ANALYSIS -eq 1 ]
    then

      ## center protein in box      
      if [ $ANALYSIS_CENT -eq 1 ]
      then
        echo 1 0 | gmx trjconv -s md.tpr -f md.xtc -o md_cent.xtc -pbc mol -center -ur compact -quiet
      fi
      
      ## calculate min and max distances from protein (reference) to membrane (selection)
      if [ $ANALYSIS_DIST -eq 1 ]
      then
        # global distances
        gmx pairdist -f md.xtc -n index.ndx -tu ns -type min -ref 1 -sel LIP -o dist_min.xvg -quiet
        gmx pairdist -f md.xtc -n index.ndx -tu ns -type max -ref 1 -sel LIP -o dist_max.xvg -quiet
      
        # residue distances (between each residue and closest lipid, format: column1: time, column2:res1-lip-dist, columnN+1=resN-lip-dist)
        gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel LIP -o dist_res_min.xvg -quiet
        gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type max -ref 1 -refgrouping res -sel LIP -o dist_res_max.xvg -quiet
	if [ $j -eq 1 ]
	then 
	  # calc the distances to specific lipids
	  gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PS -o dist_res_min_PS.xvg -quiet
	  gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PIP2 -o dist_res_min_PIP2.xvg -quiet
          gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PIP3 -o dist_res_min_PIP3.xvg -quiet
	elif [ $j -eq 2 ]
	then 
	  gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PS -o dist_res_min_PS.xvg -quiet
	  gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PIP2 -o dist_res_min_PIP2.xvg -quiet
        elif [ $j -eq 3 ]
	then 
	  gmx pairdist -f md.xtc -s md.tpr -n index.ndx -tu ns -type min -ref 1 -refgrouping res -sel PS -o dist_res_min_PS.xvg -quiet
        fi
      fi
      
      ## calculate COM distance from protein to membrane
      if [ $ANALYSIS_DIST_COM -eq 1 ]
      then
	gmx distance -f md_cent.xtc -s md.tpr -n index.ndx -select "com of group LIP plus com of group Protein" -tu ns -oall dist_com_cent -quiet
      fi

      ## calculate 2D density maps in xy plane 
      if [ $ANALYSIS_DENS -eq 1 ]
      then
        echo 1  | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PROT -quiet
        echo 13 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PC -quiet
        if [ $j -eq 1 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
          echo 15 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP2 -quiet
	  echo 16 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP3 -quiet
        fi
        if [ $j -eq 2 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
	  echo 15 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PIP2 -quiet
        fi
        if [ $j -eq 3 ]
        then
          echo 14 | gmx densmap -f md.xtc -n index.ndx -b $time -o densmap_PS -quiet
        fi
      fi
      
      ## make reference structure (tpr file) for rotmat 
      if [ $key_frame -eq 0 ] 
      then
        #using the last frame of the first rep
        if [ $j -eq $lip_for_ref -a $i -eq 0 ] 
        then
          gmx grompp -f $prod -c md.gro -r md.gro -p topol.top -o md_ref.tpr -quiet -n index.ndx
          mv md_ref.tpr ../.
        fi
      else  
        # using the frame with closest bound residue (key frame)
        # generate key frame index file
        if [ $j -eq $lip_for_ref -a $i -eq $key_frame_rep ] 
        then  
          cat << EOF > frame_index.ndx
[ frames ]

$key_frame

EOF
          # extract key frame from trajectory
          echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -quiet
	  
	  # generate tpr file from key frame
          gmx grompp -f $prod -c key_frame.gro -r key_frame.gro -p topol.top -o md_ref.tpr -quiet -n index.ndx

	  # center, for vizualization
	  echo 1 0 | gmx trjconv -s md_ref.tpr -f key_frame.gro -o key_frame_cent.gro -pbc mol -center -ur compact -quiet
          
          # convert to pdb (for visualization)
          echo 0 | gmx trjconv -f key_frame_cent.gro -s md.tpr -o key_frame_cent.pdb -quiet

	  # copy tpr file to parent directory
          mv md_ref.tpr ../.

        fi
      fi

      ## calculate rotation matrix of all frames with respect to reference structure (../md_ref.tpr)
      if [ $ANALYSIS_RZZ -eq 1 ]
      then
	echo 1 | gmx trjconv -f md.xtc -s md.tpr -o md_prot.xtc -quiet
	echo 1 | gmx trjconv -f md.gro -s md.tpr -o md_prot.gro -quiet # for visualization in vmd
        echo 1 | gmx rotmat -s ../md_ref.tpr -f md_prot.xtc -fitxy -o Rzz.xvg -quiet
      fi

      ## calculate RMSD
      if [ $ANALYSIS_RMSD -eq 1 ]
      then
        BB=$(($sel + 1))
        #echo $BB $BB | gmx rms -f eq.xtc -s eq.tpr -tu ns -quiet -o rmsd_eq.xvg -n index.ndx
        echo $BB $BB | gmx rms -f md.xtc -s md.tpr -tu ns -quiet -o rmsd_md.xvg -n index.ndx
      fi

      ## calculate solvent accessible surface area
      if [ $ANALYSIS_SASA -eq 1 ]
      then
	if [ $j -eq 2 ]
	then
	  surface=27
	  output=15
	elif [ $j -eq 5 ]
	then
	  surface=29
	  output=16
        fi
	cat << EOF > input.sasa
$surface
$output
EOF
        gmx sasa -f md.xtc -s md.tpr -n index.ndx -surface -output -probe 0.4 -o sasa.xvg -quiet < input.sasa 
        rm input.sasa
      fi
      
      ## calculate contacts between PIP2 and SOL
      if [ $ANALYSIS_CONTACT -eq 1 ]
      then
        gmx select -f md.xtc -s md.tpr -n index.ndx -select 'group SOL and within 0.5 of group PIP2' -os contact.xvg -quiet 
      fi

      ## mean square displacement (MSD)
      if [ $ANALYSIS_MSD -eq 1 ]
      then
	if [ $j -eq 2 ]
	then
          PIP2=30
        elif [ $j -eq 5 ]
        then
	  PIP2=32
        fi
        echo $PIP2 | gmx msd -f md.xtc -s md.tpr -n index.ndx -mol -tu ns -lateral z -quiet 
      fi
      
      ## make pdb of system
      echo 0 | gmx trjconv -f md.gro -s md.tpr -o md.pdb -quiet

    # end ANALYSIS if statement
    fi

    ############ FINISH #####################################################################

    ## clean up
    rm \#*

    ## navigate back to protein working directory
    cd ..

  #end loop over replicas
  done

#end loop over lipid compositions
done


## navigate back to parent directory
cd ..
