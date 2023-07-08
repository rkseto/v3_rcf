#!/bin/tcsh
if ( "$1" == "" ) then      # parentheses not strictly needed in this simple case
    echo "variable is empty"
    exit
else 
    echo "variable contains $1"
endif
hadd $1\.root ./job*/jam1ana*.root




