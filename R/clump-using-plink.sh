#!/usr/local/bin/bash

output1File=$1
RPfile=$2
measure=$3 # abf or pv

output2File=$1"-"$measure


if [ "$measure" == "abf" ]; then

/software/gapi/pkg/plink/1.0.7/bin/plink --noweb \
--bfile $RPfile \
--clump $output1File \
--clump-field neg.abf \
--clump-p1 0 \
--clump-p2 0 \
--clump-r2 0.10 \
--clump-kb 500 \
--out $output2File

gzip $output2File".clumped"

fi


if [ "$measure" == "pv" ]; then

/software/gapi/pkg/plink/1.0.7/bin/plink --noweb \
--bfile $RPfile \
--clump $output1File \
--clump-field P.value \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.10 \
--clump-kb 500 \
--out $output2File

gzip $output2File".clumped"

fi


