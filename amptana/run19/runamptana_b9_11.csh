#!/bin/tcsh
set bmin=9
set bmax=11
set setnumber=3
echo "Energy is 19 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber
mkdir -p ./jobs/'job_19GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
cd ./jobs/'job_19GeV_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
ln -sf   ../../../../../ampt/generate/19/jobs/ampt_19GeV_all_$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../../amptana3.C+\(1,-1,0,0,19,$bmin,$bmax\)
mv -f amptana3.root  amptana3_19GeV_b$bmin-$bmax\_set$setnumber.root 


