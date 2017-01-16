#/bin/bash
echo $DIR

cd /home/dailey.110/analysis_oindree
source /home/dailey.110/.bash_profile
valgrind_mem ./loop10sample 12 ~/tmp 3 3 0 5