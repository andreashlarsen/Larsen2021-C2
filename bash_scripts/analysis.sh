alchem=/sansom/s157/bioc1642/Desktop/Scripts/FEP_analysis/alchemical-analysis/alchemical_analysis/alchemical_analysis.py
py2=/sansom/s157/bioc1642/anaconda2/bin/python2.7 
dir=ANALYSIS

AXON=1
axon_dir=FEP_SMURF2_1000ns

# loop over lipids
for lip in 1 2 3 4
do
# loop over repeats
for rep in 1 2 3 4
do
  # copy from Axon
  if [ $AXON -eq 1 ]
  then
    scp -r axon:/home/bioc1642/$axon_dir/*_lip${lip}_rep${rep}_*.xvg $axon_dir
  fi

  mkdir -p $dir/lip${lip}_rep${rep}/
  for lambda in `seq 0 2 40`
  do
	if [ $AXON -eq 1 ]
	then
		cp $axon_dir/md_lip${lip}_rep${rep}_lambda${lambda}.xvg ${dir}/lip${lip}_rep${rep}
	else
		cp FEP_lip${lip}/MD_${lambda}/md_${lambda}_${rep}.xvg ${dir}/lip${lip}_rep${rep}
	fi
  done

  ls -ltr ${dir}/lip${lip}_rep${rep}
  cd ${dir}/lip${lip}_rep${rep}
  $py2 $alchem -d . -o . -q ".xvg" -p "md_lip${lip}_rep${rep}_lambda" -i 1000 -t 323 -s 5000
  cd ../..
# end loop over repeats
done
# end loop over lipids
done

