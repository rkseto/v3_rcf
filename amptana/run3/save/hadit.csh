#!/bin/tcsh
set bmin=4
set bmax=5
echo "hadd: Energy is 3 GeV bmin=" $bmin " bmax=" $bmax
hadd amptana3_3GeV_b$bmin-$bmax\_all.root  ./jobs/job*/amptana3_3GeV_b$bmin-$bmax.root



