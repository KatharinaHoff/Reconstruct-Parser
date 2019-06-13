#!/bin/bash

OUT=test.txt
FDIR=./ # adapt path to subdirectory with files
FILE_STEMS=(example/cell1.) # add more in_stems space separated

if [ -f $OUT ]; then
    echo "Output file $OUT already exists, deleting that file..."
    rm $OUT
fi

for i in "${FILE_STEMS[@]}";
do
	parse_coordinates.pl --in_stem=$FDIR$i --out=$OUT # adapt if you need to set cfg file
done

