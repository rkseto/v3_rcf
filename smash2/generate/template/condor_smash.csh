#!/bin/tcsh
setup 64b
setup ROOT 6.20.08
mkdir ./jobs/'job'$1
cd ./jobs/'job'$1
cp -d -r ../../template/* .
ln -sf ../../config.yaml config.yaml
./smash
cd ./data/0
set FILE = Particles.root
if ( -e  $FILE )  then
    echo "$FILE exists."
    mv Particles.root smashdone.root
else 
    echo "$FILE does not exist."
#    cd ../../..
#   mv -f 'job'$1 ../junk
endif
