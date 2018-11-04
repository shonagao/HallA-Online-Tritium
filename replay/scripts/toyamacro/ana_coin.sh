#!/bin/sh -f

for runlist in 00 01 02 03 04 
#for runlist in 111180_111184 111185_111190 111191_111195 111196_111200
do
  echo "$runlist"
  ./bin/ana_cointime -f runlist/Lambda${runlist} -w pdf/coin${runlist}.root -p param/f1_tuned.param >>log/${runlist}.log &
  #./bin/ana_cointime -f runlist/${runlist}.txt -w pdf/coin${runlist}.root -p param/f1_tuned.param >>log/${runlist}.log &
  #./bin/ana_cointime -f runlist/${runlist}.txt -w pdf/coin${runlist}.root -p param/tmp.param >>log/${runlist}.log &

  sleep 1
done

  sleep 20m

for runlist in 05 06 07 08 09 10
do
  echo "$runlist"
  ./bin/ana_cointime -f runlist/Lambda${runlist} -w pdf/coin${runlist}.root -p param/f1_tuned.param >>log/${runlist}.log &
  sleep 1
done
