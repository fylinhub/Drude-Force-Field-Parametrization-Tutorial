#!/bin/bash

# The x  is your staring distance*10
rm interE.dat
x=15
for i in {0..35}
do
myvar=`expr $x/10 | bc -l `     # The distance*10 should be devided here (x/10)

grep -n "Interaction" geom_$i.out |  awk -v var="$myvar" '{print var, $2 }'  | tail -n 1 >> interE.dat
x=$[$x+1]

done





