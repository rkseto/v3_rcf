#!/bin/tcsh
set bmin=0
set bmax=14
set setnumber=1
set pot="hard2SymPauli"
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
echo "input to doit1b is  file#=" $1 " centmin(file)=" 0 " centmax(file)=" 14 " potname=" $pot  " centmin="  0 " centmax=" 14 
ln -sf /star/u/rseto/work/v3/smash2/generate/3_b$bmin\-$bmax\_$pot\/jobs/smash2_3_b$bmin\-$bmax\_$pot\_all_$1.root smash2data.root
ls -l smash2data.root
root -q -b ../../../smash2ana1.C+\(6,-1,0,0,3,0,14\)
mv -f smash2ana1.root smash2_b$bmin\-$bmax\_$pot\_$1.root
exit

