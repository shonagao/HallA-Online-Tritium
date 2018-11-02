#!/bin/sh -f

for seg in 8
#for seg in 0 1 2 3 4 5 6 7 9 10 11 12 13 14 15
do
  echo "$seg"
  ./bin/mk_twc_tree -L -s ${seg} -f runlist/111167.txt -w rootfiles/LS2_${seg}.root -n 30000
  #./bin/mk_twc_tree -R -s ${seg} -f runlist/111167.txt -w rootfiles/RS2_${seg}.root -n 30000
done
