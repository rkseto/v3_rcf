#!/bin/tcsh
if ( "$5" == "" ) then      # parentheses not strictly needed in this simple case
    echo "variable is empty"
   echo " variables are e bmin bmax pot set_number nevents"
    exit
else 
    echo "variable contains $6"
endif

set e=$1
set bmin=$2
set bmax=$3
set pot=$4
set set=$5

cd job_$e\GeV_b$bmin\_$bmax\_$pot\_set$set
echo "hadd: Energy is" $e "GeV bmin=" $bmin " bmax=" $bmax "pot="$pot" set=" $set
hadd smashana3_$e\GeV_b$bmin-$bmax\_$pot\_set$set\.root  ./job*/smashana3_$e\_1M_b$bmin-$bmax\_$pot\_set$set\.root
mv -f smashana3_$e\GeV_b$bmin-$bmax\_$pot\_set$set\.root  ~/tran/smashroot/




