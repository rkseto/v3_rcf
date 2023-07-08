#!/bin/tcsh
set bmin=6
set bmax=7
set setnumber=3
echo "Energy is 27 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber
mkdir -p ./jobs/'job_27GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
cd ./jobs/'job_27GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
ln -sf   ../../../../../ampt/generate/27/jobs/ampt_27GeV_all_$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../../amptana3.C+\(1,-1,0,0,27,6.5,7.5\)
mv -f amptana3.root  amptana3_27GeV_b$bmin-$bmax\_set$setnumber.root 


