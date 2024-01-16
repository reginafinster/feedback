#!/bin/bash
CURRENTPATH="$PWD"

Gamma=0.24
N=3
for tau in {5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,16.5,17.0,17.5,18.0,18.5,19.0,19.5,20.0,20.5,21.0,21.5,22.0,22.5,23.0,23.5,24.0,24.5,25.0,25.5,26.0,26.5,27.0,27.5,28.0,28.5,29.0,29.5,30.0,30.5,31.5,32.0,32.5,33.0,33.5,34.0,34.5,35.0}
do

mkdir -p plots_N{$N}/magprofile_N{$N}_Gamma$Gamma
chmod a+rw plots_N{$N}/magprofile_N{$N}_Gamma$Gamma
cd bigplane/out
cp normalized_detector.pdf norm_detector_Gamma{$Gamma}_N{$N}.pdf
cp check.pdf check_Gamma{$Gamma}_N{$N}.pdf
cp magprofile_tau_$tau.pdf magprofile_N{$N}_Gamma{$Gamma}_tau_$tau.pdf
mv norm_detector_Gamma{$Gamma}_N{$N}.pdf ../../plots_N{$N}
mv check_Gamma{$Gamma}_N{$N}.pdf ../../plots_N{$N}
mv magprofile_N{$N}_Gamma{$Gamma}_tau_$tau.pdf ../../plots_N{$N}/magprofile_N{$N}_Gamma$Gamma
cd $CURRENTPATH

done

cp make_selected_tau_lines.py plots_N{$N}
cp make_tau_lines.py plots_N{$N}
cd bigplane/out/summary_stst
cp steadylist.dat ../../../plots_N{$N}
cd $CURRENTPATH/plots_N{$N}
python3 make_tau_lines.py
python3 make_selected_tau_lines.py


