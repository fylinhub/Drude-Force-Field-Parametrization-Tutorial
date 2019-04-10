#!/bin/bash


rm interE.dat
x=15
for i in {0..35}
do
myvar=`expr $x/10 | bc -l `

grep -n "Interaction" geom_$i.out |  awk -v var="$myvar" '{print var, $2 }'  | tail -n 1 >> interE.dat
x=$[$x+1]

done





