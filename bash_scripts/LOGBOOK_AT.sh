# load gmx etc. 
source ~/.bashrc

# user input
protein_folder=15
lip=PCPS # PC or PCPS or PCPSPIP2
rep=1

# performance
pinoffset=0

# go to CG2AT folder
cd /sansom/s157/bioc1642/Desktop/prj_C2/${protein_folder}/${lip}_${rep}/CG2AT*/FINAL

# input pdb (output from cg2at)
pdb=final_cg2at_aligned.pdb

#mdp files
mdp_folder=/sansom/s157/bioc1642/Desktop/prj_Notch/AT
min=$mdp_folder/min.mdp
nvt=$mdp_folder/nvt.mdp
npt=$mdp_folder/npt.mdp
prod=$mdp_folder/prod.mdp

# make input file
if [[ "$lip" = "PCPSPIP2" ]]
then
	cat << EOF > input
r POPC POPS POPI
name 26 LIP
q
EOF
elif [[ "$lip" = "PCPS" ]]
then    
        cat << EOF > input
r POPC POPS
name 24 LIP
q
EOF
elif [[ "$lip" = "PC" ]]
then
	cat << EOF > input
r POPC
name 22 LIP
q
EOF
else
	echo " "
	echo "ERROR: no valid lipid name given"
	echo " lip = $lip (should be PC, PCPS or PCPSPIP2)"
	echo " "
fi

gmx make_ndx -f final_cg2at_aligned.pdb -quiet < input
rm input

# nvt equilibration
gmx grompp -f $nvt -c final_cg2at_aligned.pdb -r final_cg2at_aligned.pdb -p topol_final.top -n index.ndx -o nvt.tpr -maxwarn 1 -quiet
gmx mdrun -deffnm nvt -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb auto

# npt equilibration
gmx grompp -f $npt -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_final.top -n index.ndx -o npt.tpr -maxwarn 1 -quiet
gmx mdrun -deffnm npt -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb auto

# prepare and test production run
gmx grompp -f $prod -c npt.gro -r npt.gro -t npt.cpt -p topol_final.top -n index.ndx -o prod.tpr -maxwarn 1 -quiet
gmx mdrun -deffnm prod -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb auto -nsteps 100000 # test: 0.2 ns

# production run
gmx mdrun -deffnm prod -v -quiet -pin on -pinoffset 0 -ntomp 4 -ntmpi 1 -nb auto -nsteps 20000000 # 40 ns

# extend by 60,000 ps = 60 ns (and append to previous files)
gmx convert-tpr -s prod.tpr -extend 60000 -o prod_2.tpr -quiet
mv prod.tpr prod_1.tpr
mv prod_2.tpr prod.tpr

# continue production run
gmx mdrun -deffnm prod -cpi prod.cpt -append -v -quiet -pin on -pinoffset 0 -ntomp 4 -ntmpi 1 -nb auto

# center protein in box
#echo 1 0 | gmx trjconv -s prod.tpr -f prod.xtc -o prod_cent.xtc -pbc mol -center -ur compact -quiet

# move to axon
#axon_dir=AT_${protein_folder}_${lip}_${rep}
#mkdir -p $axon_dir
#cp prod.tpr $axon_dir
#scp -r $axon_dir axon:/home/bioc1642/

# back to parent directory
cd /sansom/s157/bioc1642/Desktop/prj_C2
