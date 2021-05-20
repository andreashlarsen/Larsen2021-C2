source ~/.bashrc

RZZ=0
DIST=0
CONVERT=1

mdp_folder=/sansom/s157/bioc1642/Desktop/prj_Notch/AT
prod=$mdp_folder/prod.mdp

cat << EOF > frame_index.ndx
[ frames ]

1

EOF

if [ $CONVERT -eq 1 ]
then
	folder_frames=final_frames
	mkdir -p $folder_frames
fi

for folder in 22 27 28 18 23 25 15
#for folder in 18
do
	for lip in PC PCPS PCPSPIP2
	do
		folder_name=$folder
		if [ $folder -eq 22 ]
		then
			protein=KIBRA
			if [ "$lip" = "PC" ]
			then
				rep=17
			elif [ "$lip" = "PCPS" ]
			then
				rep=15
			elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=11
			fi
		fi
		if [ $folder -eq 27 ]
                then
			protein=PI3K_front
                        if [ "$lip" = "PC" ]
                        then
                                rep=2
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=10
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=7
                        fi
                fi
		if [ $folder -eq 28 ]
                then
			protein=PI3K_back
                        if [ "$lip" = "PC" ]
                        then
                                rep=22
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=6
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=19
                        fi
                	folder_name=27
		fi
		if [ $folder -eq 18 ]
                then
			protein=RIM2
                        if [ "$lip" = "PC" ]
                        then
                                rep=10
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=23
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=11
                        fi
                fi
		if [ $folder -eq 23 ]
                then
			protein=PTEN
                        if [ "$lip" = "PC" ]
                        then
                                rep=3
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=19
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=1
                        fi
                fi
		if [ $folder -eq 25 ]
                then
			protein=SHIP2
                        if [ "$lip" = "PC" ]
                        then
                                rep=15
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=22
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=23
                        fi
                fi
		if [ $folder -eq 15 ]
                then
			protein=Smurf2
                        if [ "$lip" = "PC" ]
                        then
                                rep=23
                        elif [ "$lip" = "PCPS" ]
                        then
                                rep=1
                        elif [ "$lip" = "PCPSPIP2" ]
                        then
                                rep=24
                        fi
                fi
		cd ${folder_name}/${lip}_${rep}/CG*/FINAL
		pwd
		
		if [ $RZZ -eq 1 ]
		then
		    # calculate Rzz with respect to frame 1
		    echo 1 | gmx trjconv -f prod_cent_short.xtc -s prod.tpr -o prod_prot.xtc -quiet
		    cp ../../../../frame_index.ndx .
		    echo 0 | gmx trjconv -f prod_cent_short.xtc -s prod.tpr -o key_frame.gro -fr frame_index.ndx -quiet
		    gmx grompp -f $prod -c key_frame.gro -r key_frame.gro -p topol_final.top -o key_frame.tpr -quiet -n index.ndx -maxwarn 1
	            echo 1 | gmx rotmat -s key_frame.tpr -f prod_prot.xtc -fitxy -o Rzz.xvg -quiet
	        fi

		if [ $DIST -eq 1 ]
		then
		    # calculated com distance
		    gmx distance -f prod_cent_short.xtc -s prod.tpr -n index.ndx -select "com of group LIP plus com of group Protein" -tu ns -oall dist_com_cent -quiet
                fi

                if [ $CONVERT -eq 1 ]
		then
                if [ "$lip" = "PCPSPIP2" ]
		then
		    # convert frame to pdb without water
		    cat <<EOF > input_index
1 | 26
name 27 PROTLIP
q
EOF
		    gmx make_ndx -f prod.gro -n index.ndx -o index_extra.ndx -quiet < input_index
		    rm input_index
		    # gmx make_ndx -f prod.gro -n index_extra.ndx -o index_extra.ndx -quiet
                    echo 1 27 | gmx trjconv -f prod.gro -s prod.tpr -n index_extra.ndx -pbc mol -ur compact -center -o ${protein}_final.pdb -quiet
		    cp ${protein}_final.pdb ../../../../$folder_frames
	        fi
                fi

		# return to parent directory
		cd ../../../..
	done
done

