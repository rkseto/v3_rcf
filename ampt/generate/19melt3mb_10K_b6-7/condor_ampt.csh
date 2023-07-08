#!/bin/tcsh
mkdir ./jobs/'job'$1
cd ./jobs/'job'$1
cp -d -r ../../template/* .
sh exec


