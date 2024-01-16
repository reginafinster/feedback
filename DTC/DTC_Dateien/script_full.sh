#!/bin/bash
CURRENTPATH='/home/rfinster2/phi_DTC/'

#for Gamma in {5,10,15}
#do
#for realization in {1,2,3,4,5,6,7,8,9,10}
#do
for phi in {80,82,84,86,88,90,92,94,96,98,100};do

RUNPATH="$CURRENTPATH/phi_$phi/realization_$realization"
mkdir -p $RUNPATH
    param=`echo "scale=3; $phi/100"|bc -l`

    
	sed "s/#phi/$param/g" < parameters_master > DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg
	
mv DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg $RUNPATH
cp qssemanybody $RUNPATH
cd $RUNPATH
#qsub -mem 10 -args DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg qssemanybody  
cd 	$CURRENTPATH
done
#done
#done



