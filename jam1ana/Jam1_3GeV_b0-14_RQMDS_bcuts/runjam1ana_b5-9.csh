#!/bin/tcsh
set bmin=5
set bmax=9
set setnumber="2_b5-9"
set pot="RQMDS"
echo "Energy is 3 GeV bmin=" $bmin " bmax=" $bmax "setnumber="$setnumber" pot="$pot
mkdir -p ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
cd ./jobs/'job_3GeV_b'$bmin'_'$bmax'_'$pot'_set'$setnumber/'job'$1
ln -sf /star/u/rseto/work/v3/jam1/generate/Jam1_3GeV_b0-14_RQMDS/run/Data/Jam1_3GeV_b0-14_RQMDS_all_$1.root jam1data.root
ls -l jam1data.root
root -q -b ../../../jam1ana5.C+\(5,-1,0,0,3,4.24,8.5\)
mv -f jam1ana5.root jam1_3_100K_b0-14_RQMDS_b$bmin-$bmax_$1.root
exit
