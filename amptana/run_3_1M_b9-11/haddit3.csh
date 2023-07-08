#!/bin/tcsh
if ( "$4" == "" ) then      # parentheses not strictly needed in this simple case
    echo "variable is empty"
   echo " variables are e bmin bmax set_number"
    exit
else 
    echo "variable contains $4"
endif

set e=$1
set bmin=$2
set bmax=$3
set set=$4
cd job_$e\_b$bmin\_$bmax\_set$set
echo "hadd: Energy is" $e "GeV bmin=" $bmin " bmax=" $bmax "set=" $set
hadd amptana3_$e\GeV_b$bmin-$bmax\_set$set\_all.root  ./job*/amptana3_$e\_b$bmin-$bmax\_set$set\.root
mv -f amptana3_$e\GeV_b$bmin-$bmax\_set$set\_all.root  ~/tran/amptroot/




