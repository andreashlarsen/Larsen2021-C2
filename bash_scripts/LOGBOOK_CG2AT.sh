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

# ensure modules and other settings from batchrc file
source ~/.bashrc

###################### expect user input from here ############################################


# select frame to CG
protein_folder=29  # select protein 
lip=2              # select lipid (1=PC:PS:PIP2:PIP3, 2=PC:PS:PIP2, 3=PC:PS, 4=PC)
rep=0              # select repetition
frame=5730         # select frame (index with 1)

# activate/deactivate (1/0) modules of script
PREPARE=0         # generate index file
BM=1              # backmapping

# Number of CPUs (for parallelization), 0 for all available
NCPUs=2

########################## to here ###################################################

# define paths to mdp files
prod=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp

# set parameters depending on protein group
NEUTRAL_PROTEIN=0 # by default assume that protein has some net charge
if [ $protein_folder -eq 8 ]
then
  protein_name=PI3K
  pdb_name=6btz
elif [ $protein_folder -eq 9 ]
then
  protein_name=FV 
  pdb_name=1czv
elif [ $protein_folder -eq 11 ]
then
  protein_name=RIM2
  pdb_name=2bwq
elif [ $protein_folder -eq 12 ]
then
  protein_name=KIBRA
  pdb_name=6fjd
elif [ $protein_folder -eq 13 ]
then
  protein_name=MFGE8
  pdb_name=2l9l
elif [ $protein_folder -eq 14 ]
then
  protein_name=SHIP2
  pdb_name=5okm
elif [ $protein_folder -eq 15 ]
then
  protein_name=SMURF2
  pdb_name=2jqz
elif [ $protein_folder -eq 16 ]
then
  protein_name=PTEN 
  pdb_name=1d5r
elif [ $protein_folder -eq 17 ]
then
  protein_name=FVIII 
  pdb_name=3hny
elif [ $protein_folder -eq 18 ]
then
  protein_name=RIM2
  pdb_name=2bwq_t
elif [ $protein_folder -eq 19 ]
then
  protein_name=DLL1
  pdb_name=4xbm
elif [ $protein_folder -eq 20 ]
then
  protein_name=DLL4
  pdb_name=5mvx
  NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 21 ]
then
  protein_name=DLL4
  pdb_name=5mvx_t
  NEUTRAL_PROTEIN=1
elif [ $protein_folder -eq 22 ]
then
  protein_name=KIBRA 
  pdb_name=6fjd_t
elif [ $protein_folder -eq 23 ]
then
  protein_name=PTEN 
  pdb_name=1d5r_t
elif [ $protein_folder -eq 24 ]
then
  protein_name=SHIP2
  pdb_name=5okm_t
elif [ $protein_folder -eq 25 ]
then
  protein_name=SHIP2
  pdb_name=5okm_tt
elif [ $protein_folder -eq 26 ]
then
  protein_name=DLL1
  pdb_name=4xbm_t
elif [ $protein_folder -eq 27 ]
then
  protein_name=PI3K
  pdb_name=6bu0
elif [ $protein_folder -eq 28 ]
then
  protein_name=SHIP2_FL
  pdb_name=5okm_FL
elif [ $protein_folder -eq 29 ]
then
  protein_name=PTEN_FL
  pdb_name=1d5r_FL
elif [ $protein_folder -eq 99 ] # for testing purpose
then
  protein_name=TEST # protein name 
  pdb_name=test # # pdb name 
  NEUTRAL_PROTEIN=1
fi

# define path to pdb file
protein=/sansom/s157/bioc1642/Desktop/prj_C2/Structures/${pdb_name}_clean.pdb # truncated

# go to protein working directory
cd $protein_folder

# loop over lipid compositions (in this case only 1)
for j in $(seq $lip $lip)
do

  # folder name
  if [ $j -eq 1 ]
  then
    Folderprefix=PCPSPIP2PIP3
  elif [ $j -eq 2 ]
  then
    Folderprefix=PCPSPIP2
  elif [ $j -eq 3 ]
  then
    Folderprefix=PCPS
  elif [ $j -eq 4 ]
  then
    Folderprefix=PC
  fi

  # loop over replicas (in this case only 1)
  for i in $(seq $rep $rep)
  do
    
    # define folder name
    folder=${Folderprefix}_$i
    
    # go to working directory
    cd $folder
    
    ############ PREPARE BACKMAPING ###########################################################
    if [ $PREPARE -eq 1 ]
    then

      # make index file
      if [ $j -eq 1 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 28 SOL
r POPC POPS POP2 POP3
name 29 LIP
a BB
r POPS
name 31 PS
r POP2
name 32 PIP2
r POP3
name 33 PIP3
q
EOF
      elif [ $j -eq 2 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 26 SOL
r POPC POPS POP2
name 27 LIP
a BB
r POPS
name 29 PS
r POP2
name 30 PIP2
q
EOF
      elif [ $j -eq 3 ]
      then 
        cat << EOF > index.input
r W WF NA+ Ion
name 24 SOL
r POPC POPS
name 25 LIP
a BB
r POPS
name 27 PS
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 1 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 16 SOL
r POPC
name 17 LIP
a BB
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 0 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 22 SOL
r POPC
name 23 LIP
a BB
q
EOF
      fi
      gmx make_ndx -f ${protein_name}_Bilayer.gro -quiet < index.input
      rm index.input

      # see overview of index file 
      # gmx make_ndx -f ${protein_name}_Bilayer.gro -n index.ndx -quiet
  
    # end PREPARE if statement
    fi

    ############ BACKMAPPING ##################################################################
    if [ $BM -eq 1 ]
    then    

      # generate key frame index file
      cat << EOF > frame_index.ndx
[ frames ]

$frame

EOF
      # extract and center key frame from trajectory
      #echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -quiet
      echo 1 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -center -pbc mol -quiet
	  
      # generate tpr file from key frame
      gmx grompp -f $prod -c key_frame.gro -r key_frame.gro -p topol.top -o key_frame.tpr -quiet -n index.ndx

      # backmapping
      if [ $NCPUs -eq 0 ]
      then
        python ~/Desktop/Scripts/cg2at_owen/cg2at/cg2at.py -c key_frame.tpr -a $protein -ff charmm36-mar2019-updated -fg martini_2-2_charmm36 -w tip3p -loc CG2AT
      else
        python ~/Desktop/Scripts/cg2at_owen/cg2at/cg2at.py -c key_frame.tpr -a $protein -ff charmm36-mar2019-updated -fg martini_2-2_charmm36 -w tip3p -ncpus $NCPUs -loc CG2AT
      fi
    fi
    
    ############ FINISH #####################################################################

    # clean up
    rm \#*

    # navigate back to protein working directory
    cd ..

  # end loop over replicas
  done

# end loop over lipid compositions
done

# navigate back to parent directory
cd ..




