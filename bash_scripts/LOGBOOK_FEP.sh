#!/bin/bash
source ~/.bashrc

#user input
protein_folder=15 # use 28 for PI3K_back
sim_time=1000 # time per frame in ns

# loop over reps
for rep in 4
do

# parallelisation with mdrun
pinoffset=0
noCPUs=2
noGPUs=1
nb=auto

# get protein name and protein_rep
if [ $protein_folder -eq 22 ]
then
	protein_name=KIBRA
	protein_rep=11
elif [ $protein_folder -eq 27 ]
then
	protein_name=PI3K
	protein_rep=7
elif [ $protein_folder -eq 28 ]
then
	protein_name=PI3K_back
	protein_rep=19
	protein_folder=27	
elif [ $protein_folder -eq 18 ]
then
        protein_name=RIM2
	protein_rep=11
elif [ $protein_folder -eq 23 ]
then
        protein_name=PTEN
	protein_rep=1
elif [ $protein_folder -eq 25 ]
then
        protein_name=SHIP2
	protein_rep=23
elif [ $protein_folder -eq 15 ]
then
        protein_name=SMURF2
	protein_rep=24
fi

# get nsteps from sim_time (with a time step of 30 fs)
nsteps=$((sim_time*1000*1000/30))

# make and go to working directory
working_dir=${protein_name}_${sim_time}ns
mkdir -p $working_dir
cd $working_dir

# make axon directory
axon_dir=FEP_${protein_name}_${sim_time}ns
mkdir -p $axon_dir

# copy (and rename) trajectory files to run example
dir_import=/sansom/s157/bioc1642/Desktop/prj_C2/$protein_folder/PCPSPIP2_$protein_rep
cp $dir_import/key_frame.gro input.gro
cp $dir_import/topol.top .
cp $dir_import/index.ndx sys.ndx
cp $dir_import/Protein.itp .

# import script files
cp /sansom/s157/bioc1642/Desktop/Scripts/FEP/FEP_examples_robin/em_RFEP.mdp .
cp /sansom/s157/bioc1642/Desktop/Scripts/FEP/FEP_examples_robin/md_RFEP.mdp .

# edit script files
sed -i -e 's/tc-grps              = PROTEIN LIPID SOL_ION/tc-grps              = Protein LIP SOL/g' md_RFEP.mdp
sed -i '1 a define           = -DFLEXIBLE' em_RFEP.mdp
sed -i '1 a define           = -DFLEXIBLE' md_RFEP.mdp
sed -i -e 's/dt                   = 0.02/dt                   = 0.03/g' md_RFEP.mdp
sed -i -e "s/nsteps               = 12500000/nsteps               = $nsteps/g" md_RFEP.mdp

# edit topology
sed -i '3 a #include "/sansom/s157/bioc1642/Desktop/Scripts/FEP/FEP_examples_robin/POP2_to_POPC_version2.itp"' topol.top
sed -i -e 's/martini_v2.2.itp/martini_v2.2_dum.itp/g' topol.top
echo "NAP              5" >> topol.top
decrese_NA=/sansom/s157/bioc1642/Desktop/prj_C2/FEP/decrese_no_of_NA_by_x.py
python $decrese_NA 5
rm topol.top
mv topol_NA.top topol.top
sed -i -e 's/POP2             4/POP2             4/g' topol.top
sed -i -e 's/POP2             4/POP2             1\nPOP2             1\nPOP2             1\nPOP2             1/' topol.top

# loop over lipid number
for lip in 1 2 3 4
do
	echo "lip = $lip"
        cp topol.top topol_lip${lip}.top
        line=$(( 20 + lip))
        sed -i "$line s/POP2             1/P2PC             1/" topol_lip${lip}.top
	
	# make a new directory to keep things neat
	dir=FEP_lip$lip
	mkdir -p $dir

	# we do lambda 0 to 40 with steps of 2. Makes it easier to "fill in" windows later on
	for i in $(seq 0 2 40)
	do
		
		if [ $rep -eq 0 ]
		then
        		########################
        		# minimise your window #
        		########################
        		mkdir -p $dir/EM_$i                                             # make a directory for minimising lambda window $i
        		sed "s/##INIT##/$i/g" em_RFEP.mdp > $dir/EM_$i/em_$i.mdp                # this sets the lambda number of your mdp 
        		gmx grompp -f $dir/EM_$i/em_$i.mdp -c input.gro -r input.gro -p topol_lip$lip.top -n sys.ndx -o $dir/EM_$i/em_$i.tpr -maxwarn 2 -quiet
        		# note that we need maxwarn here, but you should CHECK THE WARNINGS anyway. "Some parameters for bonded interaction.." is fine to ignore.
        		gmx mdrun -v -deffnm $dir/EM_$i/em_$i -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb cpu -quiet
		fi

        	###################
        	# production runs #
        	###################
        	mkdir -p $dir/MD_$i
        	sed "s/##INIT##/$i/g" md_RFEP.mdp > $dir/MD_$i/md_$i.mdp
        	for rep in $rep
        	do	
                	gmx grompp -f $dir/MD_$i/md_$i.mdp -c $dir/EM_$i/em_$i.gro -r $dir/EM_$i/em_$i.gro -p topol_lip$lip.top -n sys.ndx -o $dir/MD_$i/md_${i}_${rep}.tpr -maxwarn 2 -quiet
                	# you can either mdrun here (possibly just 50 ns per window) or scp straight over to axon etc. 12500000 steps gives 250 ns, feel free to change
                	cp $dir/MD_$i/md_${i}_${rep}.tpr $axon_dir/md_lip${lip}_rep${rep}_lambda${i}.tpr

                	#gmx mdrun -v -deffnm $dir/MD_$i/md_${i}_${rep} -nsteps 1666667 -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb -quiet 

        	# end loop over rep
        	done

	# end loop over lambda/i
	done

# end loop over lips	
done

# copy content to axon
if [ $rep -eq 0 ]
then
	scp -r $axon_dir axon:/home/bioc1642
else
	scp -r $axon_dir/*_rep${rep}_*.tpr axon:/home/bioc1642/$axon_dir
fi

#end loop over rep
done
