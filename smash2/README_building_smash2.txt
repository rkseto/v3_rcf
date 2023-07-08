updated for smash2  July 19, 2022
On pythia build - kill -march=native

using eigen, gsl, from old smash1 build. new pythia though.

rm CMakeCache.txt
 cmake .. -DPythia_CONFIG_EXECUTABLE=/star/u/rseto/work/v3/smash2/pythia8307/bin/pythia8-config -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_PREFIX_PATH=/star/u/rseto/work/v3/smash/smashtobnl/eigen/ -DGSL_ROOT_DIR=/star/u/rseto/work/v3/smash/smashtobnl/gsl -DUSE_HEPMC=OFF

maybe keep root so you can do root output
before cmake

setup 64b
setup ROOT 6.20.08 
(from Xiaoyu Liu - mattermost)

this must be done to run smash -so  make sure it is in all .csh file

Not sure this is kosher but in smash/3rdparty/virtest/vir/test.h I commented out a like with 
print(xs)      and
(extra_data)...
for extra_data comment out noinline, global, b. Traits, })

saved as the original version as test.h.sav

in all CMakeLists.txt kill =native  (MARCH)
then do the cmake and build
Then check CMakeOutput.log and CMakeError.log  and make sure "native" is never used
---------------------------------------
instructions for smash1

ln -s ~/work/v3/smash/smashtobnl/smash-master_copy/* .
mkdir build
cd build
 cmake .. -DPythia_CONFIG_EXECUTABLE=/star/u/rseto/work/v3/smash/smashtobnl/pythia8303/bin/pythia8-config -DCMAKE_PREFIX_PATH=/star/u/rseto/work/v3/smash/smashtobnl/eigen/ -DGSL_ROOT_DIR=/star/u/rseto/work/v3/smash/smashtobnl/gsl -DUSE_ROOT=OFF -DUSE_HEPMC=OFF
make
smash



rm CMakeCache.txt
 cmake .. -DPythia_CONFIG_EXECUTABLE=/star/u/rseto/work/v3/smash/smashtobnl/pythia8303/bin/pythia8-config -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_PREFIX_PATH=/star/u/rseto/work/v3/smash/smashtobnl/eigen/ -DGSL_ROOT_DIR=/star/u/rseto/work/v3/smash/smashtobnl/gsl -DUSE_ROOT=OFF -DUSE_HEPMC=OFF


in CMakeCache.txt
Kill CXX_CMARCH=Native
