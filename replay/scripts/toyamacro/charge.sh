#!/bin/sh -f

if [ $# -ne 2 ] ; then
  echo "Usage : ./charge.sh [startrunnum] [endrunnum]"
  echo "eg : ./charge.sh 111400 111500"
  exit 1
fi

rm -f scaler_charge.txt

irun=$1
frun=$2

run=$irun

while [ $run -le $frun ]
  do

echo $run
root -l -b << EOF
.L charge.C
charge($run)
EOF

  sleep 1
run=`echo "$run+1" | bc`

done

sleep 1
echo "Total charge was calculated form run$irun to run$frun"
mv scaler_charge.txt charge_info/scaler_charge_${irun}_${frun}.txt
echo "charge_info/scaler_charge_${irun}_${frun}.txt was saved"
