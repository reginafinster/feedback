#!/bin/bash
CURRENTPATH="$PWD"


for realization in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
do
for phi in {98,100}
do
for cutoff in {4,5,6,7,8,9,10,11}
do

RUNPATH="$CURRENTPATH/cutoff_1e-$cutoff/phi_$phi/realization_$realization"
mkdir -p $RUNPATH
param=`echo "scale=3; $phi/100"|bc -l`
param2=`echo "scale=3; $cutoff"|bc -l`
    
sed "s/#phi/$param/g" < parameters_master > parameters_master_temp.cfg
sed "s/#cutoff/$param2/g" < parameters_master_temp.cfg > DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg

rm parameters_master_temp.cfg
mv DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg $RUNPATH
cp qssemanybody $RUNPATH
cd $RUNPATH
qsub -mem 6 -args DTC_shared_noperturb_open_feedback_phi${param}_random_int.cfg qssemanybody  
sleep 1
cd 	$CURRENTPATH
done
done
done



