#!bin/bash

psswrd=$1
mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONrunlist > MARATHONrunlist_`date +"%m-%d-%Y"`.sql


mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONanalysis > MARATHONanalysis_`date +"%m-%d-%Y"`.sql


mysqldump -u triton-user -p${psswrd} -h halladb triton-work TargetInfo > TargetInfo_`date +"%m-%d-%Y"`.sql

mysqldump -u triton-user -p${psswrd} -h halladb triton-work coda > coda_`date +"%m-%d-%Y"`.sql

mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONrunlist > MARATHONrunlist_`date +"%m-%d-%Y"`.txt


mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONanalysis > MARATHONanalysis_`date +"%m-%d-%Y"`.txt


mysqldump -u triton-user -p${psswrd} -h halladb triton-work TargetInfo > TargetInfo_`date +"%m-%d-%Y"`.txt

mysqldump -u triton-user -p${psswrd} -h halladb triton-work coda > coda_`date +"%m-%d-%Y"`.txt

#echo 'SELECT * from MARATHONrunlist' | mysql -B -u triton-user -p${psswrd} halladb/triton-work/

#echo 'SELECT * from MARATHONanalysis' | mysql -B -u triton-user -p${psswrd} halladb/triton-work/

#echo 'SELECT * from TargetInfo' | mysql -B -u triton-user -p${psswrd} halladb/triton-work/

#echo 'SELECT * from coda' | mysql -B -u triton-user -p${psswrd} halladb/triton-work/

