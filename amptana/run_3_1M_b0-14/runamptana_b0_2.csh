#!/bin/tcsh
set bmin=0
set bmax=2
set setnumber=4
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
ln -sf   ../../../../../ampt/generate/3_1M_b0-14/jobs/ampt_3GeV_all_$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../../amptana3.C+\(1,-1,0,0,3,$bmin,$bmax\)
mv -f amptana3.root  amptana3_3GeV_b$bmin-$bmax\_set$setnumber.root 


