#!/bin/tcsh
mkdir ./jobs/'job'$1
cd ./jobs/'job'$1
cp -d -r ../../template/* .
ln -sf ../../config.yaml config.yaml
./smash
cd ./data/0
set FILE = particle_lists.oscar
if ( -e  $FILE )  then
    echo "$FILE exists."
    mv particle_lists.oscar smashdone.oscar
else 
    echo "$FILE does not exist."
    cd ../../..
    mv -f 'job'$1 ../junk
endif
