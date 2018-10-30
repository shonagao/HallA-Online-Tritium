#!/bin/sh -f

for seg in 0 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
  echo "$seg"
  ./bin/s2_twc -f rootfiles/LS2_${seg}.root -w pdf/TwcLS2_${seg}.pdf -L
done
