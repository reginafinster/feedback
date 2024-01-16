#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

for tau in {13,14,15,16,17}
do
cd tau$tau
python linecatcher_phi0.py
python linecatcher_phi0.5.py
python linecatcher_phi1.py
cd $CURRENTPATH
done
