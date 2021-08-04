#!/bin/bash

BASEDIR=/var/scratch/jb/Geoffrey_rice_blast/Data
for i in 1- 11 12 12- 13 14 15 15- 16 16- 17 18 21- 22 22- 23- 24 25 27 29 32 32- 34 36 38 39- 4- 41 43 47 49 5 5- 52 54 55 56 58 59 6 60 61 64 65 69 7 7- 70 71 9
do
  echo -e "\n\nIsolate ${i}...\n"
  cd ${BASEDIR}/Isolate${i}
done

