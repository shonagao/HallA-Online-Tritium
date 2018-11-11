#!/bin/sh -f

for seg in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
  echo "$seg"
#Fbus
#  ./bin/s2_twc -f rootfiles/LS2_${seg}.root -w pdf/TwcLS2_${seg}.pdf -L

#F1TDC
./bin/s2f1_twc -f pdf/f1_test.root -w pdf/TwcLS2F1_${seg}.pdf -L -s ${seg} -n 100000
done
