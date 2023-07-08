#!/bin/bash
g++ flow_CQ.cpp -o flow_CQ `root-config --cflags --glibs` -I${STAR_LIB}/../obj/StRoot/StEpdUtil/ -L${STAR_LIB} -lStEpdUtil
