#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

for shifts in {0.0,2.0,3.0,5.0,10.0}
do
cd shifts$shifts/out
cp detector.pdf shifts{$shifts}_shiftsdown_t40.pdf
cp normalized_detector.pdf shifts{$shifts}_shiftsdown_t40_normalized.pdf
mv shifts{$shifts}_shiftsdown_t40.pdf $CURRENTPATH
mv shifts{$shifts}_shiftsdown_t40_normalized.pdf $CURRENTPATH
cd $CURRENTPATH
done
