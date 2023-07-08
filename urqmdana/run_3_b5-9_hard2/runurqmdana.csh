#!/bin/tcsh
set bmin=5
set bmax=9
set setnumber=3
set pot=hard2
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
ln -sf   /star/u/rseto/work/v3/urqmd/generate/3_100K_b$bmin\-$bmax\_$pot\/jobs/urqmd_3_100K_b$bmin\-$bmax\_$pot\GeV_all_$1.f14   urqmddata.f14
ls -l urqmddata.f14
root -b -q ../../../../urqmdana3.C+\(2,-1,0,0,3,$bmin,$bmax\)
mv -f urqmdana3.root  urqmdana3_3_100K_b$bmin-$bmax\_$pot\_set$setnumber.root 


