#!/bin/bash

# ./CNP sample_CNPs_1.lammpstrj 1.0 1.95 12 12 0 > sample_CNPs_1.lammpstrj_shell-1.0Ang.dat
# ./CNP sample_CNPs_1.lammpstrj 2.0 1.95 12 12 0 > sample_CNPs_1.lammpstrj_shell-2.0Ang.dat
# ./CNP sample_CNPs_1.lammpstrj 3.0 1.95 12 12 0 > sample_CNPs_1.lammpstrj_shell-3.0Ang.dat

# ./CNP sample_CNPs_2.lammpstrj 1.0 1.95 12 12 0 > sample_CNPs_2.lammpstrj_shell-1.0Ang.dat
# ./CNP sample_CNPs_2.lammpstrj 2.0 1.95 12 12 0 > sample_CNPs_2.lammpstrj_shell-2.0Ang.dat
# ./CNP sample_CNPs_2.lammpstrj 3.0 1.95 12 12 0 > sample_CNPs_2.lammpstrj_shell-3.0Ang.dat

# ./CNP sample_CNPs_3.lammpstrj 1.0 1.95 12 12 0 > sample_CNPs_3.lammpstrj_shell-1.0Ang.dat
# ./CNP sample_CNPs_3.lammpstrj 2.0 1.95 12 12 0 > sample_CNPs_3.lammpstrj_shell-2.0Ang.dat
# ./CNP sample_CNPs_3.lammpstrj 3.0 1.95 12 12 0 > sample_CNPs_3.lammpstrj_shell-3.0Ang.dat

# ./CNP sample_CNPs_4.lammpstrj 1.0 1.95 12 12 0 > sample_CNPs_4.lammpstrj_shell-1.0Ang.dat
# ./CNP sample_CNPs_4.lammpstrj 2.0 1.95 12 12 0 > sample_CNPs_4.lammpstrj_shell-2.0Ang.dat
# ./CNP sample_CNPs_4.lammpstrj 3.0 1.95 12 12 0 > sample_CNPs_4.lammpstrj_shell-3.0Ang.dat

# ./CNP sample_CNPs_5.lammpstrj 1.0 1.95 12 12 0 > sample_CNPs_5.lammpstrj_shell-1.0Ang.dat
# ./CNP sample_CNPs_5.lammpstrj 2.0 1.95 12 12 0 > sample_CNPs_5.lammpstrj_shell-2.0Ang.dat
# ./CNP sample_CNPs_5.lammpstrj 3.0 1.95 12 12 0 > sample_CNPs_5.lammpstrj_shell-3.0Ang.dat

gnuplot -e "datafile='sample_CNPs_1.lammpstrj_shell-1.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_1.lammpstrj_shell-2.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_1.lammpstrj_shell-3.0Ang.dat'" sample_CNP.gp

gnuplot -e "datafile='sample_CNPs_2.lammpstrj_shell-1.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_2.lammpstrj_shell-2.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_2.lammpstrj_shell-3.0Ang.dat'" sample_CNP.gp

gnuplot -e "datafile='sample_CNPs_3.lammpstrj_shell-1.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_3.lammpstrj_shell-2.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_3.lammpstrj_shell-3.0Ang.dat'" sample_CNP.gp

gnuplot -e "datafile='sample_CNPs_4.lammpstrj_shell-1.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_4.lammpstrj_shell-2.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_4.lammpstrj_shell-3.0Ang.dat'" sample_CNP.gp

gnuplot -e "datafile='sample_CNPs_5.lammpstrj_shell-1.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_5.lammpstrj_shell-2.0Ang.dat'" sample_CNP.gp
gnuplot -e "datafile='sample_CNPs_5.lammpstrj_shell-3.0Ang.dat'" sample_CNP.gp