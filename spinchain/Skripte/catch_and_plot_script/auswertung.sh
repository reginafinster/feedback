#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

for shifts in {0.0,2.0,3.0,5.0,10.0}
do
cp linecatcher.py $CURRENTPATH/shifts$shifts/out
cd shifts$shifts/out
python3 linecatcher.py
cd $CURRENTPATH
done
