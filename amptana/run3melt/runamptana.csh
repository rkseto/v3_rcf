#!/bin/tcsh
mkdir -p ./jobs/'job'$1
cd ./jobs/'job'$1
set bmin=4
set bmax=5
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax
ln -sf   ../../../../ampt/generate/3melt/jobs/ampt$1.dat   amptdata.dat
ls -l amptdata.dat
root -b -q ../../../amptana3.C+\(1,-1,0,0,3,$bmin,$bmax\)
mv -f amptana3.root  amptana3_3GeV_b$bmin-$bmax\_melt.root 


