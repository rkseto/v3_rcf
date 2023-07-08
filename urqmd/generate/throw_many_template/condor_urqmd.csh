#!/bin/tcsh
mkdir ./jobs/'job'$1
cd ./jobs/'job'$1
cp -d -r -f ../../template/* .
bash ./runqmdrnd.bash
mv test.f14 urqmddone.f14


