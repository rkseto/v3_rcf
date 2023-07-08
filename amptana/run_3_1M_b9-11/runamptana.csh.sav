#!/bin/tcsh
set bmin=9
set bmax=11
set setnumber=5
echo "Energy is 3melt GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber
mkdir -p ./jobs/'job_3GeVmelt_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeVmelt_b'$bmin'_'$bmax'_set'$setnumber/'job'$1
ln -sf   ../../../../../ampt/generate/3melt_1M_b9-11/jobs/ampt_3melt_b9-11GeV_all_$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../../amptana3.C+\(1,-1,0,0,3,$bmin,$bmax\)
mv -f amptana3.root  amptana3_3GeVmelt_b$bmin-$bmax\_set$setnumber.root 


