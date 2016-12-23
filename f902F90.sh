#!/bin/bash

for file in *.f90
do

  ln -s $file `basename $file .f90`.F90

  echo $file

done

