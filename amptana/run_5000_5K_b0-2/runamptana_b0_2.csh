#!/bin/tcsh
set bmin=0
set bmax=2
set setnumber=4
echo "Energy is 5000 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber
mkdir -p ./jobs/'job_5000GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
cd ./jobs/'job_5000GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
ln -sf   ../../../../../ampt/generate/5000_5K_b0-2/jobs/ampt_5000_b0-2GeV_all_$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../../amptana3.C+\(1,-1,0,0,5000,$bmin,$bmax\)
mv -f amptana3.root  amptana3_5000GeV_b$bmin-$bmax\_set$setnumber.root 


