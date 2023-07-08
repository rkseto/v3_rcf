#!/bin/tcsh
set bmin=0
set bmax=14
set setnumber=5
set pot=hard2
set numevt=1M
echo "Energy is 9 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_9GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_9GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
ln -sf   /star/u/rseto/work/v3/urqmd/generate/9_$numevt\_b$bmin\-$bmax\_$pot\/jobs/urqmd_9_$numevt\_b$bmin\-$bmax\_$pot\GeV_all_$1.f14   urqmddata.f14
ls -l urqmddata.f14
root -b -q ../../../../urqmdana3.C+\(2,-1,0,0,9.2,$bmin,$bmax\)
mv -f urqmdana3.root  urqmdana3_9_$numevt\_b$bmin-$bmax\_$pot\_set$setnumber.root 


