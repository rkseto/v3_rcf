#!/bin/tcsh
set bmin=0
set bmax=14
set setnumber=9_ipt3
set pot="RQMDS"
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
echo "input to doit1b is  file#=" $1 " centmin(file)=" 0 " centmax(file)=" 14 " potname=" RQMDS  " centmin="  0 " centmax=" 14 
ln -sf /star/u/rseto/work/v3/jam1/generate/Jam1_3GeV_b0-14_RQMDS/run/Data/Jam1_3GeV_b0-14_RQMDS_all_$1.root jam1data.root
ls -l jam1data.root
root -q -b ../../../jam1ana5.C+\(5,-1,0,0,3,0,14\)
mv -f jam1ana5.root jam1_3_100K_b0-14_RQMDS_$1.root
exit

set 4 is the one we care about

set 7 now has the vector stuff
