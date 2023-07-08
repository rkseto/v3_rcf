#!/bin/bash
g++ sophisticated_QC4_f2.C -o QC4_f2 `root-config --cflags --glibs` -I${STAR_LIB}/../obj/StRoot/StEpdUtil/ -L${STAR_LIB} -lStEpdUtil
