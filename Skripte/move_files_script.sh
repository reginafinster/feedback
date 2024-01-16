#!/bin/bash

CURRENTPATH='/home/regina/Dokumente/Masterarbeit/Code/code_aktueller_Stand/'
RESULTSPATH="/home/regina/Dokumente/Masterarbeit/Code/code_aktueller_Stand/calc_results"

for realization in {1,2}
do
for mu in {10,100}
do

param=`echo "scale=3; $mu/100"|bc -l`
FETCHPATH="$CURRENTPATH/mu_$mu/realization_$realization/results"
cd $FETCHPATH
cat steadystate.dat > steadystate_${param}_realization_${realization}.dat
mkdir -p $RESULTSPATH/mu_$mu
mv steadystate_${param}_realization_${realization}.dat $RESULTSPATH/mu_$mu

done
done

