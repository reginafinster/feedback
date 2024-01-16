#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

for Gamma in {0.0,0.1,0.24,0.62}
do
for position in {1,2,3}
do
for tau in {5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0}
do

mkdir -p plots/initstate$position/magprofile_Gamma$Gamma
cd Gamma$Gamma/initstate$position/out
cp normalized_detector.pdf norm_detector_Gamma{$Gamma}_init{$position}.pdf
cp check.pdf check_Gamma{$Gamma}_init{$position}.pdf
cp magprofile_tau_$tau.pdf magprofile_Gamma{$Gamma}_init{$position}_tau_$tau.pdf
mv norm_detector_Gamma{$Gamma}_init{$position}.pdf ../../../plots/initstate$position
mv check_Gamma{$Gamma}_init{$position}.pdf ../../../plots/initstate$position
mv magprofile_Gamma{$Gamma}_init{$position}_tau_$tau.pdf ../../../plots/initstate$position/magprofile_Gamma$Gamma
cd $CURRENTPATH

done
done
done

