#datafile = "sample_CNPs_1.lammpstrj_shell-3.0Ang.dat"

###### Run this script with the following command
# gnuplot -e "datafile='datafile_name'" sample_CNP.gp


# Read overall average values from the header by using system calls.
# Note: we use awk to split the line at ": " and take the second field.
OverallRadius = real(system("grep 'Radius of Gyration' " . datafile . " | awk -F': ' '{print $2}'"))
OverallDensity  = real(system("grep 'Density \\[g/cm³\\]:' " . datafile . " | awk -F': ' '{print $2}'"))
OverallCoord    = real(system("grep 'Coordination Number:' " . datafile . " | awk -F': ' '{print $2}'"))
OverallSp3      = real(system("grep 'sp3:' " . datafile . " | awk -F': ' '{print $2}'"))
OverallSp2      = real(system("grep 'sp2:' " . datafile . " | awk -F': ' '{print $2}'"))
OverallTemp     = real(system("grep 'Temperature \\[K\\]:' " . datafile . " | awk -F': ' '{print $2}'"))

set terminal pngcairo size 900,1500
set output sprintf("%s.png", datafile)

# Arrange the plots in a grid: 5 rows x 2 columns
set multiplot layout 5,2 title sprintf("Radial Property Profiles of Carbon Nanoparticle: Rg = %.2f nm", OverallRadius) font ",14"

#########################################################################
# Plot 1: Density [g/cm³] vs Shell Center [nm]
#set title "Mass Density"
set xlabel "Distance [nm]"
set ylabel "Density [g/cm³]"
set yrange [0:5]
set format x "%.1f"
set format y "%.2f"
set label 1 sprintf(" Avg: %.2f", OverallDensity) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 1:3 w lp pt 7 ps 0.5 lw 1.5 lt rgb "blue" title ""
unset label 1

#########################################################################
# Plot 6: Density [g/cm³] vs Shell Center [R/Rg]
#set title "Mass Density"
set xlabel "R/Rg"
set ylabel "Density [g/cm³]"
set yrange [0:5]
set label 1 sprintf(" Avg: %.2f", OverallDensity) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 2:3 w lp pt 7 ps 0.5 lw 1.5 lt rgb "blue" title ""
unset label 1

#########################################################################
# Plot 2: Coordination Number vs Shell Center [nm]
#set title "Coordination Number"
set xlabel "Distance [nm]"
set ylabel "Coordination Number"
set yrange [0:4]
set label 1 sprintf(" Avg: %.2f", OverallCoord) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 1:4 w lp pt 7 ps 0.5 lw 1.5 lt rgb "red" title ""
unset label 1

#########################################################################
# Plot 7: Coordination Number vs Shell Center [R/Rg]
#set title "Coordination Number"
set xlabel "R/Rg"
set ylabel "Coordination Number"
set yrange [0:4]
set label 1 sprintf(" Avg: %.2f", OverallCoord) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 2:4 w lp pt 7 ps 0.5 lw 1.5 lt rgb "red" title ""
unset label 1

#########################################################################
# Plot 3: sp3 vs Shell Center [nm]
#set title "sp3"
set xlabel "Distance [nm]"
set ylabel "sp3"
set yrange [0:1]
set label 1 sprintf(" Avg: %.3f", OverallSp3) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 1:5 w lp pt 7 ps 0.5 lw 1.5 lt rgb "green" title ""
unset label 1

#########################################################################
# Plot 8: sp3 vs Shell Center [R/Rg]
#set title "sp3"
set xlabel "R/Rg"
set ylabel "sp3"
set yrange [0:1]
set label 1 sprintf(" Avg: %.3f", OverallSp3) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 2:5 w lp pt 7 ps 0.5 lw 1.5 lt rgb "green" title ""
unset label 1

#########################################################################
# Plot 4: sp2 vs Shell Center [nm]
#set title "sp2"
set xlabel "Distance [nm]"
set ylabel "sp2"
set yrange [0:1]
set label 1 sprintf(" Avg: %.3f", OverallSp2) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 1:6 w lp pt 7 ps 0.5 lw 1.5 lt rgb "magenta" title ""
unset label 1

#########################################################################
# Plot 9: sp2 vs Shell Center [R/Rg]
#set title "sp2"
set xlabel "R/Rg"
set ylabel "sp2"
set yrange [0:1]
set label 1 sprintf(" Avg: %.3f", OverallSp2) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 2:6 w lp pt 7 ps 0.5 lw 1.5 lt rgb "magenta" title ""
unset label 1

#########################################################################
# Plot 5: Temperature [K] vs Shell Center [nm]
#set title "Temperature"
set xlabel "Distance [nm]"
set ylabel "Temperature [K]"
set format y "%g"
set yrange [1000:7000]
set label 1 sprintf(" Avg: %.1f K", OverallTemp) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 1:7 w lp pt 7 ps 0.5 lw 1.5 lt rgb "black" title ""
unset label 1

#########################################################################
# Plot 10: Temperature [K] vs Shell Center [R/Rg]
#set title "Temperature"
set xlabel "R/Rg"
set ylabel "Temperature [K]"
set format y "%g"
set yrange [1000:7000]
set label 1 sprintf(" Avg: %.1f K", OverallTemp) at graph 0.05, graph 0.90 tc rgb "black"
p datafile u 2:7 w lp pt 7 ps 0.5 lw 1.5 lt rgb "black" title ""
unset label 1

unset multiplot