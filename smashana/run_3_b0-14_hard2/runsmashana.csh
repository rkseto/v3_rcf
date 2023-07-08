#!/bin/tcsh
set bmin=0
set bmax=14
set setnumber=5
set pot=hard2
set numevt=1M
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
ln -sf   /star/u/rseto/work/v3/smash/generate/3_$numevt\_b$bmin\-$bmax\_$pot\/jobs/smash_3_$numevt\_b$bmin\-$bmax\_$pot\GeV_all_$1.oscar   smashdata.oscar
ls -l smashdata.oscar
root -b -q ../../../../smashana3.C+\(2,-1,0,0,3,$bmin,$bmax\)
mv -f smashana3.root  smashana3_3_$numevt\_b$bmin-$bmax\_$pot\_set$setnumber.root 


