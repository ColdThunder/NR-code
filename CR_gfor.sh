export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=16
rm *.F90
./f902F90.sh
make -f Makefile_gfortran clean
make -f Makefile_gfortran
./run
