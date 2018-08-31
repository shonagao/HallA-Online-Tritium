psswrd=3He3Hdata
mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONrunlist > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONrunlist_`date +"%m-%d-%Y"`.sql


mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONanalysis > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONanalysis_`date +"%m-%d-%Y"`.sql


mysqldump -u triton-user -p${psswrd} -h halladb triton-work TargetInfo > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/TargetInfo_`date +"%m-%d-%Y"`.sql

mysqldump -u triton-user -p${psswrd} -h halladb triton-work coda > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/coda_`date +"%m-%d-%Y"`.sql

mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONTargetInfo > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONTargetInfo_`date +"%m-%d-%Y"`.sql

#mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONrunlist > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONrunlist_`date +"%m-%d-%Y"`.txt


#mysqldump -u triton-user -p${psswrd} -h halladb triton-work MARATHONanalysis > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONanalysis_`date +"%m-%d-%Y"`.txt


#mysqldump -u triton-user -p${psswrd} -h halladb triton-work TargetInfo > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/TargetInfo_`date +"%m-%d-%Y"`.txt

#mysqldump -u triton-user -p${psswrd} -h halladb triton-work coda > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/coda_`date +"%m-%d-%Y"`.txt

echo 'SELECT * from MARATHONrunlist' | mysql -B -u triton-user -p${psswrd} -h halladb triton-work > /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONrunlist_`date +"%m-%d-%Y"`.tsv

echo 'SELECT * from MARATHONanalysis' | mysql -B -u triton-user -p${psswrd} -h halladb triton-work >/home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONanalysis_`date +"%m-%d-%Y"`.tsv

echo 'SELECT * from TargetInfo' | mysql -B -u triton-user -p${psswrd} -h halladb triton-work >/home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/TargetInfo_`date +"%m-%d-%Y"`.tsv

echo 'SELECT * from coda' | mysql -B -u triton-user -p${psswrd} -h halladb triton-work >/home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/coda_`date +"%m-%d-%Y"`.tsv

echo 'SELECT * from MARATHONTargetInfo' | mysql -B -u triton-user -p${psswrd} -h halladb triton-work >/home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup/MARATHONTargetInfo_`date +"%m-%d-%Y"`.tsv

dd=`date +"%m-%d-%Y-%T"`

cd /home/jbane/tritium/replay/HallA-Online-Tritium/SQL_backup
git add ./*

git commit -m "Automatic commit made for SQL backup on ${dd}"
git push https://jbane11:1234Qwer@github.com/jbane11/HallA-Online-Tritium marathon



