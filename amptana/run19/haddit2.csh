#!/bin/tcsh
set e=19
set bmin=9
set bmax=11
echo "hadd: Energy is" $e "GeV bmin=" $bmin " bmax=" $bmax
hadd amptana3_$e\GeV_b$bmin-$bmax\_all.root  ./job*/amptana3_$e\GeV_b$bmin-$bmax.root



